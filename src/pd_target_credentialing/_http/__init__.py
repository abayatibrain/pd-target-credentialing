"""Shared HTTP infrastructure for evidence-layer clients.

This subpackage is **private** (leading underscore) — callers outside the
package import from :mod:`pd_target_credentialing.evidence` and the IO
modules, which in turn use these helpers. The retry policy and cache
design are documented in ADR-0011 and ADR-0012.
"""

from __future__ import annotations

from pd_target_credentialing._http.cache import DiskCache, request_signature
from pd_target_credentialing._http.retry import HTTPRetryError, http_retry

__all__ = [
    "DiskCache",
    "HTTPRetryError",
    "http_retry",
    "request_signature",
]
