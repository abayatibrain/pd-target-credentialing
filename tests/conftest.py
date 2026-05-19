"""Pytest configuration and shared fixtures for pd_target_credentialing."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def fixtures_dir() -> Path:
    """Path to checked-in test fixtures (small files only)."""
    return Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def _set_global_seed() -> None:
    """Per §2.4: every test runs with a deterministic seed.

    NumPy / PyTorch seeds are set when the libraries are available; tests
    that don't pull those libraries (e.g., HTTP-client unit tests) are
    not penalised by their absence.
    """
    import random

    random.seed(0)

    try:
        import numpy as np
    except ImportError:
        pass
    else:
        np.random.seed(0)

    try:
        import torch
    except ImportError:
        pass
    else:
        torch.manual_seed(0)
