"""Bounded retry decorator for HTTP-layer calls.

Per ADR-0011: three attempts, exponential backoff (1s, 2s, 4s), then
:class:`HTTPRetryError` with the underlying exception attached. The
retry is intentionally bounded — an unbounded backoff would hang the
pipeline on a sustained outage, which is worse than a clean failure
that the dossier renderer can report.
"""

from __future__ import annotations

import logging
from collections.abc import Callable
from typing import TypeVar

from tenacity import (
    RetryError,
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

logger = logging.getLogger(__name__)

T = TypeVar("T")


class HTTPRetryError(RuntimeError):
    """Raised when all retry attempts fail.

    Attributes
    ----------
    last_exception
        The exception from the final attempt. Inspect for the actual
        HTTP status or socket error.
    """

    def __init__(self, message: str, last_exception: BaseException | None = None) -> None:
        super().__init__(message)
        self.last_exception = last_exception


def http_retry(
    func: Callable[..., T],
    *,
    attempts: int = 3,
    multiplier: float = 1.0,
    max_wait: float = 8.0,
    retry_on: tuple[type[BaseException], ...] = (Exception,),
) -> Callable[..., T]:
    """Decorate ``func`` with bounded exponential-backoff retries.

    Parameters
    ----------
    func
        The function to wrap.
    attempts
        Maximum number of attempts (including the first). Default 3.
    multiplier
        Backoff multiplier in seconds (so default delays are
        1s, 2s, 4s).
    max_wait
        Cap on per-attempt wait. Default 8 seconds.
    retry_on
        Tuple of exception types to retry on. Default ``(Exception,)``.

    Returns
    -------
    callable
        Wrapped function. On exhaustion raises :class:`HTTPRetryError`.

    Notes
    -----
    Designed for HTTP clients; not for arbitrary I/O. The default
    ``retry_on=(Exception,)`` is intentionally broad because httpx
    surfaces a wide range of exception types (HTTPStatusError,
    ReadTimeout, ConnectError, ...) and narrowing here would force
    every caller to repeat the same list.
    """
    wrapped = retry(
        stop=stop_after_attempt(attempts),
        wait=wait_exponential(multiplier=multiplier, max=max_wait),
        retry=retry_if_exception_type(retry_on),
        reraise=False,
    )(func)

    def call(*args: object, **kwargs: object) -> T:
        try:
            return wrapped(*args, **kwargs)
        except RetryError as e:
            last_exc = e.last_attempt.exception() if e.last_attempt else None
            logger.error(
                "HTTP call exhausted after %d attempts: %s",
                attempts,
                last_exc,
            )
            raise HTTPRetryError(
                f"call failed after {attempts} attempts",
                last_exception=last_exc,
            ) from e

    call.__name__ = getattr(func, "__name__", "call")
    call.__doc__ = func.__doc__
    return call
