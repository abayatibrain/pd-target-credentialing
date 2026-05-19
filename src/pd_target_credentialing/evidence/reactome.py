"""Reactome ContentService client.

Design is documented in ADR-0012. Reuses the shared `_http/` cache and
retry helpers (same as the OpenTargets client) so the failure modes are
symmetric across evidence sources.

PD-relevant pathways callers commonly want:

- **R-HSA-5205685** — Mitophagy (broad).
- **R-HSA-5205647** — Mitophagy: PINK1/Parkin-mediated (canonical).

Pathway participants come back from Reactome as gene symbols. Those
symbols are routed through :class:`HGNCResolver` (ADR-0010) so the
caller always receives approved HGNC symbols; anything that fails to
resolve is dropped with a ``WARNING``.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import httpx
from pydantic import BaseModel, Field

from pd_target_credentialing._http import (
    DiskCache,
    HTTPRetryError,
    http_retry,
    request_signature,
)
from pd_target_credentialing.io.hgnc import HGNCResolver, MatchType

logger = logging.getLogger(__name__)


REACTOME_DEFAULT_URL = "https://reactome.org/ContentService"
"""Production ContentService base URL."""

# Convenience constants for the two PD-relevant pathway IDs.
PATHWAY_MITOPHAGY = "R-HSA-5205685"
PATHWAY_PINK1_PARKIN_MITOPHAGY = "R-HSA-5205647"


class ReactomeError(RuntimeError):
    """Any error talking to Reactome that the client could not recover from."""


class Pathway(BaseModel):
    """Minimal pathway record."""

    id: str
    name: str
    species: str | None = None
    is_in_disease: bool = Field(default=False, alias="isInDisease")

    model_config = {"populate_by_name": True}


class ReactomeClient:
    """Versioned, cached, retried client for the Reactome ContentService.

    Parameters
    ----------
    cache_dir
        Directory for the on-disk response cache.
    hgnc_resolver
        An :class:`HGNCResolver`. Pathway participant symbols are routed
        through it before being returned.
    base_url
        Base URL. Defaults to the public production endpoint.
    transport
        Optional ``httpx`` transport for test injection.
    timeout
        Request timeout in seconds. Default 15.
    """

    def __init__(
        self,
        cache_dir: Path,
        hgnc_resolver: HGNCResolver,
        *,
        base_url: str = REACTOME_DEFAULT_URL,
        transport: httpx.BaseTransport | None = None,
        timeout: float = 15.0,
    ) -> None:
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._cache = DiskCache(self._cache_dir / "http")
        self._base_url = base_url.rstrip("/")
        self._resolver = hgnc_resolver
        self._client = httpx.Client(
            base_url=self._base_url,
            headers={"Accept": "application/json"},
            transport=transport,
            timeout=timeout,
        )
        self._version: str | None = None

    # -- version ----------------------------------------------------------

    def get_version(self) -> str:
        """Return the pinned Reactome release version, fetched once if needed."""
        if self._version is None:
            data = self._get("/data/database/version")
            # ContentService returns a bare integer-as-text for this endpoint.
            self._version = data.strip() if isinstance(data, str) else str(data)
            logger.info("Reactome pinned to release %s", self._version)
        return self._version

    # -- queries ----------------------------------------------------------

    def get_pathway(self, pathway_id: str) -> Pathway:
        """Fetch a pathway record by stable identifier.

        Parameters
        ----------
        pathway_id
            Reactome stable ID (e.g., ``"R-HSA-5205647"``).

        Returns
        -------
        Pathway
            Parsed pathway record.

        Raises
        ------
        ReactomeError
            If the request fails after retries or the response is malformed.
        """
        payload = self._get_json(f"/data/query/{pathway_id}")
        if not isinstance(payload, dict):
            raise ReactomeError(f"unexpected payload shape for {pathway_id}")
        return Pathway.model_validate(
            {
                "id": str(payload.get("stId") or pathway_id),
                "name": str(payload.get("displayName") or ""),
                "species": (payload.get("speciesName")),
                "isInDisease": bool(payload.get("isInDisease", False)),
            }
        )

    def get_pathway_participants(self, pathway_id: str) -> list[str]:
        """Get the HGNC-approved gene symbols participating in a pathway.

        Parameters
        ----------
        pathway_id
            Reactome stable ID.

        Returns
        -------
        list[str]
            Sorted, deduplicated approved HGNC symbols. Symbols that
            cannot be resolved are dropped (with a WARNING).
        """
        payload = self._get_json(f"/data/pathway/{pathway_id}/containedEvents")
        # The ContentService offers multiple participant endpoints; the
        # one most reliably returning gene symbols is
        # `/data/participants/{id}/referenceEntities`. We use that here.
        ref = self._get_json(f"/data/participants/{pathway_id}/referenceEntities")
        symbols_raw: set[str] = set()
        candidates = ref if isinstance(ref, list) else []
        for entry in candidates:
            if not isinstance(entry, dict):
                continue
            sym = entry.get("geneName") or entry.get("identifier")
            if isinstance(sym, list) and sym:
                sym = sym[0]
            if isinstance(sym, str) and sym:
                symbols_raw.add(sym)
        _ = payload  # containedEvents not used in this version; reserved for future
        approved: set[str] = set()
        for raw_symbol in symbols_raw:
            result = self._resolver.resolve(raw_symbol)
            if result.match_type == MatchType.NOT_FOUND:
                logger.warning(
                    "Reactome participant %r could not be resolved to HGNC; dropped",
                    raw_symbol,
                )
                continue
            if result.match_type == MatchType.MULTI_MAPPING:
                logger.warning(
                    "Reactome participant %r resolved ambiguously %s; dropped",
                    raw_symbol,
                    result.candidates,
                )
                continue
            if result.approved_symbol:
                approved.add(result.approved_symbol)
        return sorted(approved)

    # -- internals --------------------------------------------------------

    def _get(self, path: str) -> str:
        """Issue a GET; return raw text. Used for the version endpoint."""
        url = f"{self._base_url}{path}"
        sig = request_signature("GET", url)
        cached = self._cache.get(sig)
        if cached is not None:
            return cached.decode("utf-8")

        @http_retry
        def _do() -> httpx.Response:
            r = self._client.get(path)
            r.raise_for_status()
            return r

        try:
            response = _do()
        except HTTPRetryError as e:
            raise ReactomeError(f"Reactome request failed: {e}") from e
        self._cache.put(sig, response.content)
        return response.text

    def _get_json(self, path: str) -> object:
        """Issue a GET; return parsed JSON."""
        url = f"{self._base_url}{path}"
        sig = request_signature("GET", url)
        cached = self._cache.get(sig)
        if cached is not None:
            return json.loads(cached.decode("utf-8"))

        @http_retry
        def _do() -> httpx.Response:
            r = self._client.get(path)
            r.raise_for_status()
            return r

        try:
            response = _do()
        except HTTPRetryError as e:
            raise ReactomeError(f"Reactome request failed: {e}") from e
        self._cache.put(sig, response.content)
        return response.json()

    def close(self) -> None:
        """Release the underlying HTTP client."""
        self._client.close()

    def __enter__(self) -> ReactomeClient:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()
