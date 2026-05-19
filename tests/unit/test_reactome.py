"""Tests for the Reactome ContentService client. No live HTTP."""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import pytest
import respx

from pd_target_credentialing.evidence.reactome import (
    PATHWAY_PINK1_PARKIN_MITOPHAGY,
    REACTOME_DEFAULT_URL,
    ReactomeClient,
)
from pd_target_credentialing.io.hgnc import HGNC_REST_BASE, HGNCResolver

REACTOME_FIX = Path(__file__).parent.parent / "fixtures" / "http" / "reactome"
HGNC_FIX = Path(__file__).parent.parent / "fixtures" / "http" / "hgnc"


def _fix(p: Path, name: str) -> dict | list | str:
    raw = (p / name).read_text()
    if name.endswith(".txt"):
        return raw
    return json.loads(raw)


@respx.mock
def test_get_version(tmp_path: Path) -> None:
    route = respx.get(f"{REACTOME_DEFAULT_URL}/data/database/version").mock(
        return_value=httpx.Response(200, text=str(_fix(REACTOME_FIX, "version.txt")))
    )
    with (
        HGNCResolver(cache_dir=tmp_path / "hgnc") as hgnc,
        ReactomeClient(cache_dir=tmp_path / "reactome", hgnc_resolver=hgnc) as rc,
    ):
        v = rc.get_version()
    assert route.called
    assert v == "88"


@respx.mock
def test_get_pathway(tmp_path: Path) -> None:
    respx.get(f"{REACTOME_DEFAULT_URL}/data/query/{PATHWAY_PINK1_PARKIN_MITOPHAGY}").mock(
        return_value=httpx.Response(200, json=_fix(REACTOME_FIX, "pathway_5205647.json"))
    )
    with (
        HGNCResolver(cache_dir=tmp_path / "hgnc") as hgnc,
        ReactomeClient(cache_dir=tmp_path / "reactome", hgnc_resolver=hgnc) as rc,
    ):
        p = rc.get_pathway(PATHWAY_PINK1_PARKIN_MITOPHAGY)
    assert p.id == PATHWAY_PINK1_PARKIN_MITOPHAGY
    assert "PINK1/Parkin" in p.name
    assert p.species == "Homo sapiens"


@respx.mock
def test_get_pathway_participants_resolves_aliases(tmp_path: Path) -> None:
    # Reactome returns: PRKN, PINK1, OPTN, PARK2 (alias).
    # Each gets resolved via HGNC. PRKN and PARK2 both → PRKN; PINK1 → PINK1;
    # OPTN → OPTN. Final set: {PRKN, PINK1, OPTN}.
    respx.get(
        f"{REACTOME_DEFAULT_URL}/data/pathway/{PATHWAY_PINK1_PARKIN_MITOPHAGY}/containedEvents"
    ).mock(
        return_value=httpx.Response(200, json=_fix(REACTOME_FIX, "contained_events_5205647.json"))
    )
    respx.get(
        f"{REACTOME_DEFAULT_URL}/data/participants/{PATHWAY_PINK1_PARKIN_MITOPHAGY}/referenceEntities"
    ).mock(return_value=httpx.Response(200, json=_fix(REACTOME_FIX, "participants_5205647.json")))

    # Stub HGNC responses for each symbol the participants endpoint returned.
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PRKN").mock(
        return_value=httpx.Response(200, json=_fix(HGNC_FIX, "prkn_approved.json"))
    )
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PARK2").mock(
        return_value=httpx.Response(200, json=_fix(HGNC_FIX, "park2_alias.json"))
    )
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PINK1").mock(
        return_value=httpx.Response(
            200,
            json={
                "response": {
                    "numFound": 1,
                    "docs": [
                        {
                            "symbol": "PINK1",
                            "hgnc_id": "HGNC:14581",
                            "ensembl_gene_id": "ENSG00000158828",
                        }
                    ],
                }
            },
        )
    )
    respx.get(f"{HGNC_REST_BASE}/search/symbol/OPTN").mock(
        return_value=httpx.Response(
            200,
            json={
                "response": {
                    "numFound": 1,
                    "docs": [
                        {
                            "symbol": "OPTN",
                            "hgnc_id": "HGNC:17142",
                            "ensembl_gene_id": "ENSG00000123240",
                        }
                    ],
                }
            },
        )
    )

    with (
        HGNCResolver(cache_dir=tmp_path / "hgnc") as hgnc,
        ReactomeClient(cache_dir=tmp_path / "reactome", hgnc_resolver=hgnc) as rc,
    ):
        symbols = rc.get_pathway_participants(PATHWAY_PINK1_PARKIN_MITOPHAGY)

    # Deduplicated approved symbols.
    assert symbols == sorted({"PRKN", "PINK1", "OPTN"})


@respx.mock
def test_reactome_request_failure_raises(tmp_path: Path) -> None:
    respx.get(f"{REACTOME_DEFAULT_URL}/data/query/R-HSA-MISSING").mock(
        return_value=httpx.Response(404, text="not found")
    )
    from pd_target_credentialing.evidence.reactome import ReactomeError

    with (
        HGNCResolver(cache_dir=tmp_path / "hgnc") as hgnc,
        ReactomeClient(cache_dir=tmp_path / "reactome", hgnc_resolver=hgnc) as rc,
        pytest.raises(ReactomeError),
    ):
        rc.get_pathway("R-HSA-MISSING")
