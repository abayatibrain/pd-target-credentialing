"""Configuration dataclasses for pd_target_credentialing.

Every threshold or magic number used in the pipeline lives here (per §2.4 of
the Cowork brief: "No magic numbers"). Each field carries a docstring
explaining *why* this value, not just what.

The defaults reflect the values committed to in ADR-0001 and subsequent ADRs.
If you want to change a default, write a superseding ADR first.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class PathsConfig:
    """Filesystem layout. Override via env or CLI flags."""

    cache_dir: Path = field(
        default_factory=lambda: Path.home() / ".cache" / "pd_target_credentialing"
    )
    """Where downloaded raw data lives. Configurable so shared-filesystem
    users can override (§2.6)."""

    results_dir: Path = Path("results")
    """Where generated PNGs/SVGs are written. Each figure gets a
    `.caption.md` sidecar with dataset version + commit SHA (§2.6)."""


@dataclass(frozen=True)
class Thresholds:
    """Statistical thresholds. Defaults are defended in ADRs.

    Changing any of these without an ADR update is a defect.
    """

    fdr_alpha: float = 0.05
    """Benjamini-Hochberg FDR cutoff. See ADR-0001."""


@dataclass(frozen=True)
class AppConfig:
    """Top-level config bundle."""

    paths: PathsConfig = field(default_factory=PathsConfig)
    thresholds: Thresholds = field(default_factory=Thresholds)
    random_seed: int = 0xC0FFEE
    """Global seed used for NumPy, Python random, and PyTorch (§2.4)."""
