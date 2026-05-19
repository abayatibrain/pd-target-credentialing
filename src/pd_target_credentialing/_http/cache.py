"""Content-addressed on-disk HTTP response cache.

Per ADR-0011, the cache is content-addressed by SHA256 of the canonicalized
request signature (method + URL + sorted-keys JSON body). Cache hits return
the previously stored payload bytes; cache misses are populated by the
caller after the live request returns.

The cache is **deterministic** and **server-header-independent**, which is
the correct behavior for reproducing dossier results months after the
underlying API has rotated cache headers. Cache invalidation is the
caller's responsibility — typically by including a platform version
string in the request signature, so a release bump produces new keys.

Example
-------
>>> cache = DiskCache(Path("/tmp/cache"))
>>> sig = request_signature("POST", "https://api/x", body={"q": "PRKN"})
>>> hit = cache.get(sig)
>>> if hit is None:
...     payload = b'{"result": "..."}'
...     cache.put(sig, payload)
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class RequestSignature:
    """Canonical signature of an HTTP request, for cache keying.

    Two requests with the same signature produce the same SHA256 digest,
    regardless of dict key ordering in the body.
    """

    method: str
    url: str
    body_hash: str  # hex digest of the canonicalized body, or "" for no body

    @property
    def digest(self) -> str:
        """SHA256 hex digest of (method, url, body_hash) joined by U+001E."""
        s = f"{self.method.upper()}\x1e{self.url}\x1e{self.body_hash}"
        return hashlib.sha256(s.encode("utf-8")).hexdigest()


def request_signature(
    method: str,
    url: str,
    *,
    body: dict[str, Any] | str | bytes | None = None,
) -> RequestSignature:
    """Build a canonical :class:`RequestSignature` for a request.

    Parameters
    ----------
    method
        HTTP method (``"GET"``, ``"POST"``, ...). Case-insensitive.
    url
        Full request URL.
    body
        Optional request body. Dicts are JSON-canonicalized (sorted keys,
        no whitespace) before hashing so that semantically equal payloads
        cache the same way regardless of dict ordering. Bytes and strings
        are hashed as-is.

    Returns
    -------
    RequestSignature
        Hashable signature; call ``.digest`` for the SHA256 hex string.
    """
    if body is None:
        body_hash = ""
    else:
        if isinstance(body, dict):
            payload = json.dumps(body, sort_keys=True, separators=(",", ":")).encode("utf-8")
        elif isinstance(body, str):
            payload = body.encode("utf-8")
        elif isinstance(body, bytes):
            payload = body
        else:  # pragma: no cover — defensive
            raise TypeError(f"unsupported body type: {type(body).__name__}")
        body_hash = hashlib.sha256(payload).hexdigest()
    return RequestSignature(method=method.upper(), url=url, body_hash=body_hash)


class DiskCache:
    """A minimal content-addressed on-disk cache.

    Storage layout: ``<root>/<aa>/<rest>.bin`` where ``<aa>`` is the first
    two hex chars of the digest (avoids putting millions of files in one
    directory).

    The cache stores raw bytes; serialization is the caller's
    responsibility. The convention used by this repo is to store the raw
    response body (typically JSON bytes) and let the caller parse.
    """

    def __init__(self, root: Path) -> None:
        """Create a cache rooted at ``root``.

        Parameters
        ----------
        root
            Directory to use. Created if it does not exist.
        """
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)

    def _path_for(self, sig: RequestSignature) -> Path:
        digest = sig.digest
        return self.root / digest[:2] / f"{digest[2:]}.bin"

    def get(self, sig: RequestSignature) -> bytes | None:
        """Return the cached payload for ``sig``, or ``None`` on miss.

        Parameters
        ----------
        sig
            Request signature.

        Returns
        -------
        bytes or None
            The cached body bytes, or ``None`` if the cache does not
            contain this signature.
        """
        path = self._path_for(sig)
        if path.exists():
            logger.debug("disk cache HIT: %s", sig.digest)
            return path.read_bytes()
        logger.debug("disk cache MISS: %s", sig.digest)
        return None

    def put(self, sig: RequestSignature, payload: bytes) -> None:
        """Store ``payload`` under ``sig``.

        Parameters
        ----------
        sig
            Request signature.
        payload
            Raw bytes to store.
        """
        path = self._path_for(sig)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        logger.debug("disk cache PUT: %s (%d bytes)", sig.digest, len(payload))

    def clear(self) -> None:
        """Remove all cached entries (but keep the root directory)."""
        for shard in self.root.iterdir():
            if shard.is_dir():
                for entry in shard.iterdir():
                    entry.unlink()
                shard.rmdir()
