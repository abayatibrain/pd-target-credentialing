"""OpenTargets GraphQL client.

Design is documented in ADR-0011. Three guarantees this client provides:

1. **Platform-version pinning.** The platform release version is read
   once at construction time via the ``meta { dataVersion }`` query and
   becomes part of every cache key, so a future OT release does not
   silently change a previously cached score.
2. **Content-addressed response cache.** A re-run of the same query on
   the same release is a zero-network operation.
3. **Bounded retry.** Three attempts, exponential backoff, then a clean
   :class:`OpenTargetsError`.

This module deliberately does **not** implement the score-aggregation
logic. That is ADR-0008's territory and is pending Armin sign-off. Here
we only fetch and return the raw evidence channels.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import httpx
from pydantic import BaseModel, Field

from pd_target_credentialing._http import (
    DiskCache,
    HTTPRetryError,
    http_retry,
    request_signature,
)

logger = logging.getLogger(__name__)


OPENTARGETS_DEFAULT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
"""Production GraphQL endpoint."""


class OpenTargetsError(RuntimeError):
    """Any error talking to OpenTargets that the client could not recover from.

    Includes both HTTP-layer failures (exhausted retries) and
    GraphQL-layer errors returned in the response ``errors`` field.
    """


@dataclass(frozen=True)
class PlatformVersion:
    """The pinned platform release information.

    Attributes
    ----------
    data_version
        Free-form version string published by OT (``"24.06"`` and so on).
    api_version
        Free-form API version published by OT (kept for the audit trail).
    """

    data_version: str
    api_version: str


# ----------------------------------------------------------------------
# Pydantic response models
# ----------------------------------------------------------------------


class TargetSummary(BaseModel):
    """Minimal target record returned by the symbol lookup query."""

    id: str
    """OpenTargets Ensembl-anchored target ID (e.g., ``"ENSG00000145335"``)."""
    # Field name matches the GraphQL field; the camelCase deviation from
    # PEP 8 is deliberate and the N815 ruff rule is not enabled.
    approvedSymbol: str
    biotype: str | None = None


class DatatypeScore(BaseModel):
    """One per-datatype score for a target-disease pair."""

    id: str
    """OT datatype identifier (e.g., ``"genetic_association"``)."""
    score: float


class AssociationByDatatype(BaseModel):
    """Result of an association-by-datatype query."""

    target_id: str
    disease_id: str
    overall_score: float
    by_datatype: list[DatatypeScore] = Field(default_factory=list)


class TractabilityBucket(BaseModel):
    """One tractability assessment row."""

    modality: str
    """e.g., ``"smallMolecule"``, ``"antibody"``, ``"PROTAC"``."""
    value: bool
    label: str | None = None


# ----------------------------------------------------------------------
# Client
# ----------------------------------------------------------------------


class OpenTargetsClient:
    """Versioned, cached, retried GraphQL client for OpenTargets.

    Parameters
    ----------
    cache_dir
        Directory for the on-disk response cache (per ADR-0011).
    base_url
        GraphQL endpoint. Defaults to the public production endpoint.
    transport
        Optional ``httpx`` transport for test injection.
    timeout
        Request timeout in seconds. Default 15.
    """

    def __init__(
        self,
        cache_dir: Path,
        *,
        base_url: str = OPENTARGETS_DEFAULT_URL,
        transport: httpx.BaseTransport | None = None,
        timeout: float = 15.0,
    ) -> None:
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._cache = DiskCache(self._cache_dir / "http")
        # Note: do NOT pass base_url to httpx.Client — httpx appends a
        # trailing slash when the request path is empty, which silently
        # diverges from the URL we register in the cache signature and
        # the URL respx mocks see. Always pass the full URL explicitly.
        self._base_url = base_url
        self._client = httpx.Client(
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
            },
            transport=transport,
            timeout=timeout,
        )
        self._version: PlatformVersion | None = None

    # -- version pinning --------------------------------------------------

    def version(self) -> PlatformVersion:
        """Return the pinned platform version, fetching it once if needed."""
        if self._version is None:
            query = "{ meta { name dataVersion { year month iteration } apiVersion { x y z } } }"
            payload = self._post(query, {}, include_version_in_key=False)
            meta = (payload.get("data") or {}).get("meta") or {}
            dv = meta.get("dataVersion") or {}
            av = meta.get("apiVersion") or {}
            self._version = PlatformVersion(
                data_version=f"{dv.get('year', 0)}.{dv.get('month', 0):02d}"
                if isinstance(dv.get("month"), int)
                else str(dv.get("year") or "unknown"),
                api_version=".".join(str(av.get(k, "?")) for k in ("x", "y", "z")),
            )
            logger.info(
                "OpenTargets pinned: data=%s api=%s",
                self._version.data_version,
                self._version.api_version,
            )
        return self._version

    # -- queries ----------------------------------------------------------

    def get_target_by_symbol(self, symbol: str) -> TargetSummary | None:
        """Look up a target by approved gene symbol.

        Parameters
        ----------
        symbol
            HGNC-approved gene symbol. Callers must resolve aliases via
            :class:`pd_target_credentialing.io.hgnc.HGNCResolver` first.

        Returns
        -------
        TargetSummary or None
            ``None`` if the symbol does not resolve at OT.
        """
        query = """
        query TargetBySymbol($q: String!) {
          search(queryString: $q, entityNames: ["target"]) {
            hits {
              id
              entity
              object { ... on Target { id approvedSymbol biotype } }
            }
          }
        }
        """
        payload = self._post(query, {"q": symbol})
        hits = ((payload.get("data") or {}).get("search") or {}).get("hits") or []
        for hit in hits:
            obj = hit.get("object") or {}
            if (
                hit.get("entity") == "target"
                and obj.get("approvedSymbol", "").upper() == symbol.upper()
            ):
                return TargetSummary.model_validate(obj)
        return None

    def get_association_by_datatype(self, target_id: str, disease_id: str) -> AssociationByDatatype:
        """Get per-datatype association scores for a target-disease pair.

        Parameters
        ----------
        target_id
            OT target ID (e.g., ``"ENSG00000145335"`` for SNCA).
        disease_id
            EFO ID for the disease (e.g., ``"EFO_0002508"`` for PD).

        Returns
        -------
        AssociationByDatatype
            Overall score and per-datatype breakdown. Always returns —
            an OT response with no matching association yields zeros.
        """
        query = """
        query Association($t: String!, $d: String!) {
          target(ensemblId: $t) {
            associatedDiseases(efoIds: [$d]) {
              rows {
                disease { id }
                score
                datatypeScores { id score }
              }
            }
          }
        }
        """
        payload = self._post(query, {"t": target_id, "d": disease_id})
        rows = (
            ((payload.get("data") or {}).get("target") or {}).get("associatedDiseases") or {}
        ).get("rows") or []
        if not rows:
            return AssociationByDatatype(
                target_id=target_id,
                disease_id=disease_id,
                overall_score=0.0,
                by_datatype=[],
            )
        row = rows[0]
        return AssociationByDatatype(
            target_id=target_id,
            disease_id=disease_id,
            overall_score=float(row.get("score") or 0.0),
            by_datatype=[
                DatatypeScore.model_validate(d) for d in (row.get("datatypeScores") or [])
            ],
        )

    def get_tractability(self, target_id: str) -> list[TractabilityBucket]:
        """Tractability buckets for a target.

        Parameters
        ----------
        target_id
            OT target ID.

        Returns
        -------
        list[TractabilityBucket]
            Per-modality assessment. Empty list if OT has no tractability
            data for this target.
        """
        query = """
        query Tractability($t: String!) {
          target(ensemblId: $t) {
            tractability { id modality label value }
          }
        }
        """
        payload = self._post(query, {"t": target_id})
        raw = ((payload.get("data") or {}).get("target") or {}).get("tractability") or []
        out: list[TractabilityBucket] = []
        for item in raw:
            if not isinstance(item, dict):
                continue
            out.append(
                TractabilityBucket(
                    modality=str(item.get("modality") or item.get("id") or "unknown"),
                    value=bool(item.get("value")),
                    label=item.get("label"),
                )
            )
        return out

    # -- internals --------------------------------------------------------

    def _post(
        self,
        query: str,
        variables: dict[str, Any],
        *,
        include_version_in_key: bool = True,
    ) -> dict[str, Any]:
        """Issue a GraphQL POST with cache + retry, return parsed JSON."""
        body: dict[str, Any] = {"query": query, "variables": variables}
        if include_version_in_key:
            body["_pin"] = self.version().data_version
        sig = request_signature("POST", self._base_url, body=body)
        cached = self._cache.get(sig)
        if cached is not None:
            return json.loads(cached.decode("utf-8"))

        @http_retry
        def _do_call() -> httpx.Response:
            r = self._client.post(
                self._base_url,
                json={"query": query, "variables": variables},
            )
            r.raise_for_status()
            return r

        try:
            response = _do_call()
        except HTTPRetryError as e:
            raise OpenTargetsError(f"OT request failed: {e}") from e

        payload_bytes = response.content
        payload = json.loads(payload_bytes.decode("utf-8"))
        if payload.get("errors"):
            raise OpenTargetsError(f"OT GraphQL errors: {payload['errors']}")
        self._cache.put(sig, payload_bytes)
        return payload

    def close(self) -> None:
        """Release the underlying HTTP client."""
        self._client.close()

    def __enter__(self) -> OpenTargetsClient:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.close()
