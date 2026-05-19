"""Loader for Kamath et al. 2022 human midbrain DA atlas.

ADR-0001 ratifies Kamath 2022 as the **cross-cohort validator** for
DA-anchored claims. ADR-0005 Q5.3 puts DA subtypes (SOX6+ vulnerable
vs CALB1+ less-vulnerable) **in scope for v1.0.0**, with Kamath as the
reference atlas for the cross-check.

Mode design mirrors :mod:`pd_target_credentialing.io.smajic2022`:

- ``"toy"`` reuses the committed toy fixture for shape/schema testing.
  Subtype labels (``SOX6_pos`` / ``CALB1_pos`` / ``unresolved``) are
  *synthesized* in toy mode purely so downstream code that consumes
  subtypes has something to consume; the labels are not biologically
  meaningful and the loader logs a warning to that effect.
- ``"real"`` will download from the Single Cell Portal (Broad Institute,
  SCP1768). Not implemented in this commit; tracked for v0.2.0.

Citation
--------
Kamath T. et al. (2022) *Nat Neurosci* 25:588-595.
doi:[10.1038/s41593-022-01061-1](https://doi.org/10.1038/s41593-022-01061-1)
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    import anndata as ad

from pd_target_credentialing.io.smajic2022 import TOY_FIXTURE

logger = logging.getLogger(__name__)


DATASET_NAME = "kamath2022"
SCP_ACCESSION = "SCP1768"

DA_SUBTYPE_LABELS: tuple[str, ...] = ("SOX6_pos", "CALB1_pos", "unresolved")
"""DA-neuron subtypes used by ADR-0005 Q5.3's cross-check."""


class KamathLoadError(RuntimeError):
    """Anything that prevents the loader from returning a valid AnnData."""


def load_kamath2022(
    *,
    mode: Literal["real", "toy"] = "toy",
    cache_dir: Path | None = None,
    fixture_path: Path | None = None,
) -> ad.AnnData:
    """Load the Kamath 2022 midbrain DA atlas as a cross-cohort reference.

    Parameters
    ----------
    mode
        ``"toy"`` reads the toy fixture and synthesizes DA-subtype labels.
        ``"real"`` downloads from the Broad Single Cell Portal (not yet
        implemented).
    cache_dir
        Override for the default cache directory (real mode only).
    fixture_path
        Override for the toy fixture location.

    Returns
    -------
    AnnData
        Same schema as the Smajic loader, plus an additional ``obs``
        column ``da_subtype`` with values in :data:`DA_SUBTYPE_LABELS`.
        For non-DA cells, ``da_subtype`` is ``"unresolved"``.

    Raises
    ------
    KamathLoadError
        Toy fixture missing.
    NotImplementedError
        Real-mode download not yet implemented.
    """
    if mode == "real":
        raise NotImplementedError(
            "Real-mode Kamath 2022 loader not implemented yet. "
            f"Single Cell Portal accession: {SCP_ACCESSION}. "
            "Use mode='toy' for now."
        )
    if mode != "toy":
        raise ValueError(f"unknown mode: {mode!r}")
    return _load_toy(fixture_path)


def _load_toy(fixture_path: Path | None) -> ad.AnnData:
    """Read the toy fixture and synthesize DA-subtype labels for tests."""
    import anndata as ad

    path = fixture_path if fixture_path is not None else TOY_FIXTURE
    if not path.exists():
        raise KamathLoadError(
            f"Toy fixture not found at {path}. Generate it via "
            "`python tests/fixtures/make_toy_anndata.py`."
        )
    adata = ad.read_h5ad(path)
    # Synthesize subtype labels deterministically so tests are reproducible.
    # The labels here are NOT biologically meaningful — they are just enough
    # plumbing for the downstream subtype cross-check to be exercised.
    logger.warning(
        "Kamath toy-mode loader: synthesizing DA-subtype labels for plumbing "
        "tests. The labels are NOT biological and must not be reported."
    )
    subtypes = _synthesize_subtypes(adata)
    adata.obs["da_subtype"] = subtypes
    adata.uns["kamath_subtype_source"] = "synthesized-toy-mode"
    return adata


def _synthesize_subtypes(adata: ad.AnnData) -> list[str]:
    """Assign each cell a subtype label.

    DA neurons get split 60/40 SOX6+/CALB1+ based on a deterministic hash
    of the cell barcode. Non-DA cells get "unresolved".
    """
    out: list[str] = []
    for idx, row in adata.obs.iterrows():
        ct_field = row.get("celltype_true", row.get("celltype", "unknown"))
        if str(ct_field) == "DA_neuron":
            # Deterministic split based on cell index hash.
            digit = int(str(idx).split("_")[-1]) if "_" in str(idx) else hash(str(idx))
            out.append("SOX6_pos" if (digit % 10) < 6 else "CALB1_pos")
        else:
            out.append("unresolved")
    return out
