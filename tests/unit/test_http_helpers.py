"""Tests for the shared _http helpers (cache + retry)."""
from __future__ import annotations

from pathlib import Path

import pytest

from pd_target_credentialing._http import (
    DiskCache,
    HTTPRetryError,
    http_retry,
    request_signature,
)


# ---------- request_signature ------------------------------------------


def test_signature_is_method_case_insensitive() -> None:
    a = request_signature("GET", "https://x/y")
    b = request_signature("get", "https://x/y")
    assert a.digest == b.digest


def test_signature_canonicalizes_dict_body_key_order() -> None:
    a = request_signature("POST", "https://x", body={"a": 1, "b": 2})
    b = request_signature("POST", "https://x", body={"b": 2, "a": 1})
    assert a.digest == b.digest


def test_signature_differs_on_url_change() -> None:
    a = request_signature("GET", "https://x/a")
    b = request_signature("GET", "https://x/b")
    assert a.digest != b.digest


def test_signature_differs_on_body_change() -> None:
    a = request_signature("POST", "https://x", body={"q": "A"})
    b = request_signature("POST", "https://x", body={"q": "B"})
    assert a.digest != b.digest


def test_signature_rejects_unknown_body_type() -> None:
    with pytest.raises(TypeError):
        request_signature("POST", "https://x", body=12345)  # type: ignore[arg-type]


def test_signature_handles_string_and_bytes_bodies() -> None:
    s = request_signature("POST", "https://x", body="hello")
    b = request_signature("POST", "https://x", body=b"hello")
    assert s.digest == b.digest


# ---------- DiskCache --------------------------------------------------


def test_disk_cache_miss_then_hit(tmp_path: Path) -> None:
    cache = DiskCache(tmp_path)
    sig = request_signature("GET", "https://x/y")
    assert cache.get(sig) is None
    cache.put(sig, b"hello")
    assert cache.get(sig) == b"hello"


def test_disk_cache_clear(tmp_path: Path) -> None:
    cache = DiskCache(tmp_path)
    sig = request_signature("GET", "https://x")
    cache.put(sig, b"v")
    assert cache.get(sig) == b"v"
    cache.clear()
    assert cache.get(sig) is None


def test_disk_cache_creates_root(tmp_path: Path) -> None:
    nested = tmp_path / "a" / "b" / "c"
    DiskCache(nested)
    assert nested.exists()


# ---------- http_retry -------------------------------------------------


def test_http_retry_returns_value_on_success() -> None:
    @http_retry
    def f() -> int:
        return 42

    assert f() == 42


def test_http_retry_recovers_after_transient_failures() -> None:
    attempts = {"n": 0}

    @http_retry
    def f() -> str:
        attempts["n"] += 1
        if attempts["n"] < 3:
            raise RuntimeError("flaky")
        return "ok"

    assert f() == "ok"
    assert attempts["n"] == 3


def test_http_retry_raises_after_exhaustion() -> None:
    @http_retry
    def f() -> None:
        raise RuntimeError("always fails")

    with pytest.raises(HTTPRetryError) as exc:
        f()
    assert isinstance(exc.value.last_exception, RuntimeError)
