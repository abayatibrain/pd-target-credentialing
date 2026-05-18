"""Tests for the HGNC resolver. No live HTTP — all calls are respx-mocked."""
from __future__ import annotations

import json
import logging
from pathlib import Path

import httpx
import pytest
import respx

from pd_target_credentialing.io.hgnc import (
    HGNC_REST_BASE,
    HGNCResolver,
    MatchType,
    MultiMappingError,
)

FIXTURES = Path(__file__).parent.parent / "fixtures" / "http" / "hgnc"


def _load_fixture(name: str) -> dict[str, object]:
    return json.loads((FIXTURES / name).read_text())


# ---------- happy paths ------------------------------------------------


@respx.mock
def test_resolve_approved_symbol(tmp_path: Path) -> None:
    route = respx.get(f"{HGNC_REST_BASE}/search/symbol/PRKN").mock(
        return_value=httpx.Response(200, json=_load_fixture("prkn_approved.json"))
    )
    with HGNCResolver(cache_dir=tmp_path) as r:
        result = r.resolve("PRKN")
    assert route.called
    assert result.match_type == MatchType.APPROVED
    assert result.approved_symbol == "PRKN"
    assert result.hgnc_id == "HGNC:8607"
    assert result.ensembl_id == "ENSG00000185345"


@respx.mock
def test_resolve_alias_emits_warning(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PARK2").mock(
        return_value=httpx.Response(200, json=_load_fixture("park2_alias.json"))
    )
    with caplog.at_level(logging.WARNING, logger="pd_target_credentialing.io.hgnc"):
        with HGNCResolver(cache_dir=tmp_path) as r:
            result = r.resolve("PARK2")
    assert result.match_type == MatchType.ALIAS
    assert result.approved_symbol == "PRKN"
    # The substitution must surface at WARNING (per §3.4).
    assert any("alias substitution" in rec.message.lower() for rec in caplog.records)


@respx.mock
def test_resolve_is_case_insensitive(tmp_path: Path) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PRKN").mock(
        return_value=httpx.Response(200, json=_load_fixture("prkn_approved.json"))
    )
    with HGNCResolver(cache_dir=tmp_path) as r:
        a = r.resolve("prkn")
        b = r.resolve("  PRKN  ")
    assert a.approved_symbol == b.approved_symbol == "PRKN"


# ---------- error paths ------------------------------------------------


@respx.mock
def test_resolve_not_found(tmp_path: Path) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/NOSUCHGENE").mock(
        return_value=httpx.Response(200, json=_load_fixture("not_found.json"))
    )
    with HGNCResolver(cache_dir=tmp_path) as r:
        result = r.resolve("NOSUCHGENE")
    assert result.match_type == MatchType.NOT_FOUND
    assert result.approved_symbol is None


@respx.mock
def test_resolve_multi_mapping_does_not_silently_pick(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/AMBIG").mock(
        return_value=httpx.Response(200, json=_load_fixture("multi_mapping.json"))
    )
    with caplog.at_level(logging.ERROR, logger="pd_target_credentialing.io.hgnc"):
        with HGNCResolver(cache_dir=tmp_path) as r:
            result = r.resolve("AMBIG")
    assert result.match_type == MatchType.MULTI_MAPPING
    assert result.approved_symbol is None
    assert sorted(result.candidates) == ["GENE1", "GENE2"]
    assert any("multi-mapping" in rec.message.lower() for rec in caplog.records)


@respx.mock
def test_resolve_strict_raises_on_multi_mapping(tmp_path: Path) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/AMBIG").mock(
        return_value=httpx.Response(200, json=_load_fixture("multi_mapping.json"))
    )
    with HGNCResolver(cache_dir=tmp_path) as r:
        with pytest.raises(MultiMappingError):
            r.resolve_strict("AMBIG")


@respx.mock
def test_resolve_strict_raises_keyerror_on_not_found(tmp_path: Path) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/NOSUCHGENE").mock(
        return_value=httpx.Response(200, json=_load_fixture("not_found.json"))
    )
    with HGNCResolver(cache_dir=tmp_path) as r:
        with pytest.raises(KeyError):
            r.resolve_strict("NOSUCHGENE")


# ---------- caching ----------------------------------------------------


@respx.mock
def test_disk_cache_avoids_second_network_call(tmp_path: Path) -> None:
    route = respx.get(f"{HGNC_REST_BASE}/search/symbol/PRKN").mock(
        return_value=httpx.Response(200, json=_load_fixture("prkn_approved.json"))
    )
    # First resolver hits the network.
    with HGNCResolver(cache_dir=tmp_path) as r:
        r.resolve("PRKN")
    assert route.call_count == 1

    # A fresh resolver pointed at the same cache must not call again.
    with HGNCResolver(cache_dir=tmp_path) as r:
        result = r.resolve("PRKN")
    assert route.call_count == 1  # unchanged
    assert result.approved_symbol == "PRKN"


# ---------- multi-mapping CSV side effect -------------------------------


@respx.mock
def test_multimapping_writes_csv(tmp_path: Path) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/AMBIG").mock(
        return_value=httpx.Response(200, json=_load_fixture("multi_mapping.json"))
    )
    csv_path = tmp_path / "review.csv"
    with HGNCResolver(cache_dir=tmp_path, multimapping_csv=csv_path) as r:
        r.resolve("AMBIG")
    assert csv_path.exists()
    rows = csv_path.read_text().strip().splitlines()
    assert rows[0] == "input_symbol,candidates,resolved_at"
    # Header + one data row.
    assert rows[1].startswith("AMBIG,GENE1;GENE2,")


# ---------- batch ------------------------------------------------------


@respx.mock
def test_resolve_many_preserves_inputs(tmp_path: Path) -> None:
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PRKN").mock(
        return_value=httpx.Response(200, json=_load_fixture("prkn_approved.json"))
    )
    respx.get(f"{HGNC_REST_BASE}/search/symbol/PARK2").mock(
        return_value=httpx.Response(200, json=_load_fixture("park2_alias.json"))
    )
    with HGNCResolver(cache_dir=tmp_path) as r:
        out = r.resolve_many(["PRKN", "park2"])
    assert set(out.keys()) == {"PRKN", "PARK2"}
    assert out["PARK2"].approved_symbol == "PRKN"
    assert out["PRKN"].approved_symbol == "PRKN"
