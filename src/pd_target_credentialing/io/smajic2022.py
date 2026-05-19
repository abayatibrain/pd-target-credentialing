"""Loader for Smajic et al. 2022 substantia nigra snRNA-seq.

ADR-0001 ratifies Smajic 2022 as the primary dataset.

This module has two modes:

- **Production** (``mode="real"``): downloads from GEO (GSE178265) with
  SHA256 verification, caches under ``$XDG_CACHE_HOME/pd_target_credentialing/``.
  Multi-gigabyte; not exercised in CI.
- **Toy** (``mode="toy"``): loads the synthetic 500-nucleus fixture
  committed under ``tests/fixtures/toy_smajic.h5ad``. Returns an AnnData
  with the exact same schema as the real mode, so downstream modules
  cannot tell them apart at runtime.

Production mode is intentionally a placeholder in v0.1.0 — the download
+ SHA256 plumbing is the next implementation slice (post-foundations).
For now, real-mode raises ``NotImplementedError`` with a clear message,
and the toy-mode path is the one tests exercise.

Citation
--------
Smajic S. et al. (2022) *Brain* 145(3):964-978.
doi:[10.1093/brain/awab406](https://doi.org/10.1093/brain/awab406)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import anndata as ad

logger = logging.getLogger(__name__)


DATASET_NAME = "smajic2022"
GEO_ACCESSION = "GSE178265"
TOY_FIXTURE = Path(__file__).resolve().parents[3] / "tests" / "fixtures" / "toy_smajic.h5ad"


class SmajicLoadError(RuntimeError):
    """Anything that prevents the loader from returning a valid AnnData."""


def load_smajic2022(
    *,
    mode: Literal["real", "toy"] = "toy",
    cache_dir: Path | None = None,
    fixture_path: Path | None = None,
) -> ad.AnnData:
    """Load the Smajic 2022 substantia nigra snRNA-seq dataset.

    Parameters
    ----------
    mode
        ``"toy"`` (default) reads the committed synthetic fixture, suitable
        for tests and for the demo notebook's smoke run. ``"real"`` downloads
        the published cohort from GEO; not implemented in this commit.
    cache_dir
        Override for the default cache directory (real mode only).
    fixture_path
        Override for the toy fixture location. Only useful for tests that
        want to point at a different synthetic fixture.

    Returns
    -------
    AnnData
        Raw counts in ``.X`` (and mirrored in ``.layers["counts"]``).
        Required ``obs`` columns: ``donor_id``, ``condition``, ``age``,
        ``sex``, ``PMI``. Required ``var`` columns: ``gene_name``,
        ``is_mito``.

    Raises
    ------
    SmajicLoadError
        If toy mode is requested and the fixture is missing.
    NotImplementedError
        If real mode is requested. The download plumbing lands in a
        subsequent commit.
    """
    if mode == "real":
        # ADR-0001 + brief Section 2.6: real-data download is gated behind a
        # SHA256-verified, cache-aware loader. That implementation lands in
        # the v0.2.0 milestone (see STATUS.md "Plan for next week").
        raise NotImplementedError(
            "Real-mode Smajic 2022 loader not implemented yet. Use "
            "mode='toy' for tests and the demo notebook smoke run. "
            f"GEO accession (for the future fetcher): {GEO_ACCESSION}. "
            "Cache dir would be: "
            f"{cache_dir or '$XDG_CACHE_HOME/pd_target_credentialing/'}"
        )

    if mode != "toy":
        raise ValueError(f"unknown mode: {mode!r}")

    return _load_toy(fixture_path)


def _load_toy(fixture_path: Path | None) -> ad.AnnData:
    """Read the committed toy fixture and verify its schema."""
    import anndata as ad

    path = fixture_path if fixture_path is not None else TOY_FIXTURE
    if not path.exists():
        raise SmajicLoadError(
            f"Toy fixture not found at {path}. "
            "Generate it by running `python tests/fixtures/make_toy_anndata.py`."
        )
    adata = ad.read_h5ad(path)
    _verify_schema(adata)
    logger.info(
        "Loaded toy Smajic fixture: %d nuclei x %d genes from %s",
        adata.n_obs,
        adata.n_vars,
        path,
    )
    return adata


def _verify_schema(adata: ad.AnnData) -> None:
    """Check the loaded AnnData has the columns downstream modules expect."""
    required_obs = {"donor_id", "condition", "age", "sex", "PMI"}
    missing = required_obs - set(adata.obs.columns)
    if missing:
        raise SmajicLoadError(f"Smajic AnnData is missing required obs columns: {sorted(missing)}")
    if "is_mito" not in adata.var.columns:
        raise SmajicLoadError("Smajic AnnData is missing required var column: 'is_mito'")
    if "counts" not in adata.layers:
        # Tolerate counts-only-in-X fixtures by copying X to the counts layer
        # so downstream pseudobulk DE can always reach raw counts.
        logger.warning("Smajic AnnData has no 'counts' layer; copying .X to .layers['counts'].")
        adata.layers["counts"] = adata.X.copy()
