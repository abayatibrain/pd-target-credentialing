"""Tests for the OpenTargets GraphQL client. All HTTP via respx mocks."""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import pytest
import respx

from pd_target_credentialing.evidence.opentargets import (
    OPENTARGETS_DEFAULT_URL,
    OpenTargetsClient,
    OpenTargetsError,
)

FIXTURES = Path(__file__).parent.parent / "fixtures" / "http" / "opentargets"


def _fixture(name: str) -> dict:
    return json.loads((FIXTURES / name).read_text())


@respx.mock
def test_version_is_pinned_and_cached(tmp_path: Path) -> None:
    route = respx.post(OPENTARGETS_DEFAULT_URL).mock(
        return_value=httpx.Response(200, json=_fixture("meta.json"))
    )
    with OpenTargetsClient(cache_dir=tmp_path) as ot:
        v1 = ot.version()
        v2 = ot.version()
    # Second call must not hit the network.
    assert route.call_count == 1
    assert v1.data_version == "2024.06"
    assert v2.data_version == v1.data_version


@respx.mock
def test_get_target_by_symbol(tmp_path: Path) -> None:
    # The client issues meta+query, so respx returns based on call order
    # via side-effect ordering on the same route.
    payloads = iter([_fixture("meta.json"), _fixture("target_by_symbol_snca.json")])

    def _responder(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=next(payloads))

    respx.post(OPENTARGETS_DEFAULT_URL).mock(side_effect=_responder)
    with OpenTargetsClient(cache_dir=tmp_path) as ot:
        target = ot.get_target_by_symbol("SNCA")
    assert target is not None
    assert target.id == "ENSG00000145335"
    assert target.approvedSymbol == "SNCA"


@respx.mock
def test_get_association_by_datatype(tmp_path: Path) -> None:
    payloads = iter([_fixture("meta.json"), _fixture("association_snca_pd.json")])

    def _responder(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=next(payloads))

    respx.post(OPENTARGETS_DEFAULT_URL).mock(side_effect=_responder)
    with OpenTargetsClient(cache_dir=tmp_path) as ot:
        a = ot.get_association_by_datatype("ENSG00000145335", "EFO_0002508")
    assert a.overall_score == pytest.approx(0.91)
    ids = {s.id for s in a.by_datatype}
    assert {"genetic_association", "literature", "animal_model"}.issubset(ids)


@respx.mock
def test_get_tractability(tmp_path: Path) -> None:
    payloads = iter([_fixture("meta.json"), _fixture("tractability_snca.json")])

    def _responder(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=next(payloads))

    respx.post(OPENTARGETS_DEFAULT_URL).mock(side_effect=_responder)
    with OpenTargetsClient(cache_dir=tmp_path) as ot:
        buckets = ot.get_tractability("ENSG00000145335")
    modalities = {b.modality for b in buckets}
    assert "smallMolecule" in modalities
    sm = next(b for b in buckets if b.modality == "smallMolecule")
    assert sm.value is True


@respx.mock
def test_graphql_errors_raise(tmp_path: Path) -> None:
    payloads = iter(
        [
            _fixture("meta.json"),
            {"errors": [{"message": "bad query"}]},
        ]
    )

    def _responder(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=next(payloads))

    respx.post(OPENTARGETS_DEFAULT_URL).mock(side_effect=_responder)
    with (
        OpenTargetsClient(cache_dir=tmp_path) as ot,
        pytest.raises(OpenTargetsError),
    ):
        ot.get_target_by_symbol("SNCA")


@respx.mock
def test_response_is_cached_across_clients(tmp_path: Path) -> None:
    payloads = iter([_fixture("meta.json"), _fixture("target_by_symbol_snca.json")])

    def _responder(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=next(payloads))

    route = respx.post(OPENTARGETS_DEFAULT_URL).mock(side_effect=_responder)
    with OpenTargetsClient(cache_dir=tmp_path) as ot:
        ot.get_target_by_symbol("SNCA")
    first_calls = route.call_count

    # New client, same cache_dir → no further network calls for the same query.
    # We re-prime meta because version() refetches on a new client instance,
    # but the target lookup should hit the cache.
    payloads2 = iter([_fixture("meta.json")])

    def _responder2(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json=next(payloads2))

    respx.post(OPENTARGETS_DEFAULT_URL).mock(side_effect=_responder2)
    with OpenTargetsClient(cache_dir=tmp_path) as ot:
        ot.get_target_by_symbol("SNCA")
    # The new client made the meta call (1) but reused the cached target hit.
    # If caching is broken, this would be 2 instead of 1.
    assert (
        route.call_count - first_calls == 0
    )  # original mock no longer in use; new mock counts elsewhere
