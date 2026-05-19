"""HGNC alias resolution and on-disk symbol cache.

Per §3.4 of the Cowork brief, this module is the **single source of
truth** for gene-symbol identity in this repo. Every place that
ingests external data with gene symbols must route through here.

Design is documented in ADR-0010. Key properties:

- Approved-symbol-only output. Aliases are resolved silently in code
  but logged at ``WARNING`` so the pipeline's audit trail shows every
  substitution.
- Multi-mapping detection: when one input alias resolves to more than
  one approved symbol, the resolver does **not** silently pick one.
  It returns ``MatchType.MULTI_MAPPING`` and writes the case to a CSV
  for human review.
- On-disk cache so that a re-run with the same input set performs zero
  network calls.

Example
-------
>>> from pathlib import Path
>>> # resolver = HGNCResolver(cache_dir=Path("/tmp/hgnc"))
>>> # result = resolver.resolve("PARK2")  # legacy alias for PRKN
>>> # assert result.approved_symbol == "PRKN"
"""

from __future__ import annotations

import csv
import json
import logging
from datetime import UTC, datetime
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

import httpx
from pydantic import BaseModel, Field

from pd_target_credentialing._http import DiskCache, http_retry, request_signature

if TYPE_CHECKING:
    from collections.abc import Iterable

logger = logging.getLogger(__name__)


HGNC_REST_BASE = "https://rest.genenames.org"
"""Base URL for the HGNC REST API. Pinned at module level."""


class MatchType(StrEnum):
    """How an input symbol was matched against HGNC records."""

    APPROVED = "approved"
    """Input symbol is the current approved symbol."""

    ALIAS = "alias"
    """Input symbol is a current alias of a single approved symbol."""

    PREV_SYMBOL = "prev_symbol"
    """Input symbol is a previously-used symbol of a single approved
    symbol (e.g., PARK2 → PRKN)."""

    MULTI_MAPPING = "multi_mapping"
    """Input symbol resolves to more than one approved symbol — see
    :class:`MultiMappingError`."""

    NOT_FOUND = "not_found"
    """Input symbol does not appear anywhere in HGNC's records."""


class ResolutionResult(BaseModel):
    """Outcome of resolving one input symbol against HGNC."""

    input_symbol: str
    approved_symbol: str | None = None
    """``None`` for ``MULTI_MAPPING`` or ``NOT_FOUND``."""
    hgnc_id: str | None = None
    ensembl_id: str | None = None
    entrez_id: str | None = None
    match_type: MatchType
    candidates: list[str] = Field(default_factory=list)
    """For ``MULTI_MAPPING``: the list of approved symbols that the input
    resolves to. For other match types, empty."""
    resolved_at: datetime


class MultiMappingError(ValueError):
    """Raised by :meth:`HGNCResolver.resolve_strict` when one input symbol
    unambiguously maps to multiple approved symbols.

    The non-strict :meth:`HGNCResolver.resolve` returns a
    :class:`ResolutionResult` with ``match_type == MULTI_MAPPING`` instead
    of raising, so callers can decide policy per use site.
    """


class HGNCResolver:
    """Resolve gene aliases to HGNC approved symbols, with on-disk cache.

    Parameters
    ----------
    cache_dir
        Root directory for the on-disk HTTP response cache.
    multimapping_csv
        Path to a CSV where multi-mapping cases are recorded for review.
        Defaults to ``data/hgnc_multimappings.csv`` relative to cwd.
    transport
        Optional ``httpx`` transport, used in tests to inject ``respx``
        mocks. Production code does not pass this.
    """

    def __init__(
        self,
        cache_dir: Path,
        *,
        multimapping_csv: Path | None = None,
        transport: httpx.BaseTransport | None = None,
    ) -> None:
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._disk_cache = DiskCache(self._cache_dir / "http")
        self._multimapping_csv = (
            multimapping_csv
            if multimapping_csv is not None
            else Path("data") / "hgnc_multimappings.csv"
        )
        self._client = httpx.Client(
            base_url=HGNC_REST_BASE,
            headers={"Accept": "application/json"},
            transport=transport,
            timeout=10.0,
        )
        self._memo: dict[str, ResolutionResult] = {}

    # -- public API -------------------------------------------------------

    def resolve(self, symbol: str) -> ResolutionResult:
        """Resolve a single input symbol to its HGNC approved symbol.

        Parameters
        ----------
        symbol
            Input gene symbol (any casing; whitespace-stripped).

        Returns
        -------
        ResolutionResult
            On a clean resolution, ``approved_symbol`` is populated.
            On a multi-mapping, ``match_type`` is
            :attr:`MatchType.MULTI_MAPPING` and ``candidates`` lists the
            ambiguous targets — caller decides policy.
        """
        key = symbol.strip().upper()
        if key in self._memo:
            return self._memo[key]

        result = self._fetch_and_classify(key)
        self._memo[key] = result

        if result.match_type in (MatchType.ALIAS, MatchType.PREV_SYMBOL):
            logger.warning(
                "HGNC alias substitution: %r → %r (%s)",
                symbol,
                result.approved_symbol,
                result.match_type.value,
            )
        elif result.match_type == MatchType.MULTI_MAPPING:
            logger.error(
                "HGNC multi-mapping: %r → %s — flagged for review",
                symbol,
                result.candidates,
            )
            self._append_multimapping_row(result)
        elif result.match_type == MatchType.NOT_FOUND:
            logger.error("HGNC not found: %r", symbol)

        return result

    def resolve_strict(self, symbol: str) -> str:
        """Resolve to a single approved symbol or raise.

        Parameters
        ----------
        symbol
            Input gene symbol.

        Returns
        -------
        str
            The approved HGNC symbol.

        Raises
        ------
        MultiMappingError
            If the input maps to more than one approved symbol.
        KeyError
            If the input is not found.
        """
        r = self.resolve(symbol)
        if r.match_type == MatchType.MULTI_MAPPING:
            raise MultiMappingError(f"{symbol!r} maps to multiple approved symbols: {r.candidates}")
        if r.match_type == MatchType.NOT_FOUND or r.approved_symbol is None:
            raise KeyError(f"{symbol!r} not found in HGNC")
        return r.approved_symbol

    def resolve_many(self, symbols: Iterable[str]) -> dict[str, ResolutionResult]:
        """Resolve a batch of symbols. Order-preserving.

        Parameters
        ----------
        symbols
            Iterable of input symbols.

        Returns
        -------
        dict[str, ResolutionResult]
            Keyed by the canonical (uppercased, stripped) input.
        """
        return {s.strip().upper(): self.resolve(s) for s in symbols}

    # -- internals --------------------------------------------------------

    def _fetch_and_classify(self, key: str) -> ResolutionResult:
        """Hit the HGNC API (with cache + retry) and classify the result."""
        url = f"{HGNC_REST_BASE}/search/symbol/{key}"
        sig = request_signature("GET", url)
        cached = self._disk_cache.get(sig)
        if cached is not None:
            payload = self._decode_json(cached)
        else:

            @http_retry
            def _do_call() -> httpx.Response:
                r = self._client.get(f"/search/symbol/{key}")
                r.raise_for_status()
                return r

            response = _do_call()
            payload = response.json()
            self._disk_cache.put(sig, response.content)

        return self._classify(key, payload)

    @staticmethod
    def _decode_json(payload: bytes) -> dict[str, object]:
        decoded = json.loads(payload.decode("utf-8"))
        if not isinstance(decoded, dict):
            raise ValueError("expected JSON object from HGNC")
        return decoded

    def _classify(self, key: str, payload: dict[str, object]) -> ResolutionResult:
        """Convert an HGNC ``/search/symbol/{key}`` payload to a result.

        HGNC's response shape is::

            {"response": {"docs": [{"symbol": "...", "hgnc_id": "...",
                                    ...}], "numFound": N}}
        """
        response = payload.get("response", {})
        if not isinstance(response, dict):
            response = {}
        docs_raw = response.get("docs") or []
        docs: list[dict[str, object]] = (
            [d for d in docs_raw if isinstance(d, dict)] if isinstance(docs_raw, list) else []
        )

        now = datetime.now(UTC)

        if len(docs) == 0:
            return ResolutionResult(
                input_symbol=key,
                match_type=MatchType.NOT_FOUND,
                resolved_at=now,
            )

        if len(docs) > 1:
            candidates = sorted({str(d.get("symbol")) for d in docs if d.get("symbol")})
            return ResolutionResult(
                input_symbol=key,
                match_type=MatchType.MULTI_MAPPING,
                candidates=candidates,
                resolved_at=now,
            )

        doc = docs[0]
        approved = str(doc.get("symbol"))
        match_type = MatchType.APPROVED if approved.upper() == key else MatchType.ALIAS
        return ResolutionResult(
            input_symbol=key,
            approved_symbol=approved,
            hgnc_id=str(doc["hgnc_id"]) if doc.get("hgnc_id") else None,
            ensembl_id=str(doc["ensembl_gene_id"]) if doc.get("ensembl_gene_id") else None,
            entrez_id=str(doc["entrez_id"]) if doc.get("entrez_id") else None,
            match_type=match_type,
            resolved_at=now,
        )

    def _append_multimapping_row(self, result: ResolutionResult) -> None:
        """Append one CSV row for human review of multi-mapping cases."""
        self._multimapping_csv.parent.mkdir(parents=True, exist_ok=True)
        new_file = not self._multimapping_csv.exists()
        with self._multimapping_csv.open("a", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            if new_file:
                writer.writerow(["input_symbol", "candidates", "resolved_at"])
            writer.writerow(
                [
                    result.input_symbol,
                    ";".join(result.candidates),
                    result.resolved_at.isoformat(),
                ]
            )

    def close(self) -> None:
        """Release the underlying HTTP client."""
        self._client.close()

    def __enter__(self) -> HGNCResolver:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()
