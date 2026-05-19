"""Marker-based cell-type annotation for substantia nigra nuclei.

Implements ADR-0005's hybrid scheme (Armin sign-off 2026-05-18, Q5.1-Q5.2):

1. **Score** every nucleus against each cell type by averaging the
   log1p-normalized expression of that type's marker panel.
2. **Assign** the highest-scoring cell type.
3. **Flag ambiguous** nuclei whose top-two scores are within
   ``ambiguity_margin`` (default 0.15 per ADR-0005 Q5.2) of each other;
   those nuclei are kept in the AnnData with
   ``celltype_ambiguous=True`` and excluded from downstream DE per
   ADR-0005's consequences.

The Kamath cross-check for DA-subtype labels lives in
:mod:`pd_target_credentialing.annotate.da_subtypes`.

Example
-------
>>> # annotated = annotate_celltypes(adata)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

from pd_target_credentialing.annotate.markers import get_panel

if TYPE_CHECKING:
    import anndata as ad
    import numpy as np

    from pd_target_credentialing.annotate.markers import MarkerPanel

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AnnotationConfig:
    """Parameters for the annotation pass. Defaults from ADR-0005."""

    ambiguity_margin: float = 0.15
    """Top-two-score margin below which a nucleus is flagged ambiguous
    (ADR-0005 Q5.2). Ambiguous nuclei are excluded from downstream DE."""

    use_log1p: bool = True
    """If True, apply log1p to the count matrix before scoring (so marker
    expression contributes in normalized space). When the AnnData already
    has a log1p layer named ``log1p`` we use that instead."""

    scale_to: float = 1e4
    """Library-size scaling target before log1p (per ADR-0003)."""


def _get_log1p_matrix(adata: ad.AnnData, config: AnnotationConfig) -> np.ndarray:
    """Return a dense numpy matrix of log1p-normalized counts."""
    import numpy as np
    import scipy.sparse as sp

    if "log1p" in adata.layers:
        X = adata.layers["log1p"]
    else:
        X = adata.X
        # Ensure raw-counts shape; library-size scale then log1p
        X = X.toarray() if sp.issparse(X) else np.asarray(X)
        if config.use_log1p:
            total = X.sum(axis=1, keepdims=True)
            total[total == 0] = 1.0
            X = np.log1p(X / total * config.scale_to)
    return np.asarray(X, dtype=np.float64)


def annotate_celltypes(
    adata: ad.AnnData,
    *,
    config: AnnotationConfig | None = None,
    panel: MarkerPanel | None = None,
) -> ad.AnnData:
    """Annotate nuclei by marker-score voting.

    Parameters
    ----------
    adata
        Input AnnData (post-QC). The function returns a copy with three
        new ``obs`` columns:

        - ``celltype`` (str): assigned label, or ``"ambiguous"``.
        - ``celltype_confidence`` (float): the top score, in [0, 1] after
          per-cell rescaling so the top score is 1.0 for the most
          confident cell.
        - ``celltype_ambiguous`` (bool): True if the top-two margin is
          smaller than ``config.ambiguity_margin``.

    config
        Annotation parameters. Defaults to ADR-0005's values.
    panel
        Override the canonical panel (tests). Defaults to ``get_panel()``.

    Returns
    -------
    AnnData
        Annotated copy.
    """
    import numpy as np

    cfg = config if config is not None else AnnotationConfig()
    p = panel if panel is not None else get_panel()
    log1p = _get_log1p_matrix(adata, cfg)

    gene_to_col = {str(g): i for i, g in enumerate(adata.var.index)}
    n_cells = adata.n_obs

    # Score each (cell, cell_type) by averaging the markers present in the data.
    cell_types = p.cell_types
    scores = np.zeros((n_cells, len(cell_types)), dtype=np.float64)
    missing_markers: dict[str, list[str]] = {}
    for k, entry in enumerate(p):
        cols = [gene_to_col[m] for m in entry.markers if m in gene_to_col]
        if not cols:
            missing_markers[entry.cell_type] = list(entry.markers)
            continue
        # Mean log1p across the markers that exist
        scores[:, k] = log1p[:, cols].mean(axis=1)

    if missing_markers:
        logger.warning(
            "Cell types with NO markers present in the data (will never be assigned): %s",
            sorted(missing_markers.keys()),
        )

    # Per-cell rescale so confidence is comparable across cells.
    max_per_cell = scores.max(axis=1, keepdims=True)
    max_per_cell[max_per_cell == 0] = 1.0
    rescaled = scores / max_per_cell

    # Assign and check ambiguity
    sorted_idx = np.argsort(-rescaled, axis=1)
    top_idx = sorted_idx[:, 0]
    top_score = rescaled[np.arange(n_cells), top_idx]
    second_idx = sorted_idx[:, 1] if rescaled.shape[1] > 1 else top_idx
    second_score = rescaled[np.arange(n_cells), second_idx]

    margin = top_score - second_score
    ambiguous = margin < cfg.ambiguity_margin

    labels = np.array([cell_types[i] for i in top_idx], dtype=object)
    labels[ambiguous] = "ambiguous"

    out = adata.copy()
    out.obs["celltype"] = labels
    out.obs["celltype_confidence"] = top_score.astype(float)
    out.obs["celltype_ambiguous"] = ambiguous

    n_amb = int(ambiguous.sum())
    logger.info(
        "Annotation: %d nuclei -> %d unambiguous, %d ambiguous (margin < %.2f).",
        n_cells,
        n_cells - n_amb,
        n_amb,
        cfg.ambiguity_margin,
    )
    return out
