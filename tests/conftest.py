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
    """Per §2.4: every test runs with a deterministic seed."""
    import random

    import numpy as np

    random.seed(0)
    np.random.seed(0)
