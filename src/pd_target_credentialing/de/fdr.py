"""Multiple-testing correction for cell-type-conditioned DE.

ADR-0007 (Armin sign-off 2026-05-18) commits to:

- **Primary**: Benjamini-Hochberg FDR within each cell type's p-value
  vector independently, at alpha = 0.05.
- **Secondary**: BH-FDR across the full pooled p-value vector (all
  cell types x all genes). Reported alongside but not the headline.

Genes significant under BOTH corrections receive a "strong evidence"
flag used by the score module (ADR-0008).

Implementation is pure scipy/numpy; no statsmodels dependency.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

    from pd_target_credentialing.de.pydeseq2_runner import DEResult

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class FDRConfig:
    """FDR-correction knobs (ADR-0007-ratified defaults)."""

    alpha: float = 0.05
    """Q7.1 — pre-registered alpha."""

    pvalue_col: str = "pvalue"
    """Column in each per-celltype DE table holding raw p-values."""

    within_padj_col: str = "padj_within_celltype"
    """Output column for per-cell-type BH-FDR."""

    global_padj_col: str = "padj_global"
    """Output column for global BH-FDR (Q7.2)."""

    strong_evidence_col: str = "strong_evidence"
    """True when both within-celltype and global FDR pass alpha."""


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Vectorised Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    pvalues
        Raw p-values. NaNs are passed through.

    Returns
    -------
    np.ndarray
        BH-adjusted p-values, same shape as input.
    """
    import numpy as np

    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    if n == 0:
        return p
    nan_mask = np.isnan(p)
    if nan_mask.all():
        return p

    valid_idx = np.where(~nan_mask)[0]
    p_valid = p[valid_idx]
    order = np.argsort(p_valid)
    ranked = p_valid[order]
    m = len(ranked)
    adjusted = ranked * m / (np.arange(1, m + 1))
    # Enforce monotonicity (BH step-up procedure)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    # Unsort
    out_valid = np.empty(m, dtype=float)
    out_valid[order] = adjusted
    out = np.full(n, np.nan, dtype=float)
    out[valid_idx] = out_valid
    return out


def apply_fdr(
    de_results: dict[str, DEResult], config: FDRConfig | None = None
) -> dict[str, DEResult]:
    """Add per-celltype and global FDR columns to every DE table.

    Parameters
    ----------
    de_results
        Output of :func:`run_pydeseq2_per_celltype`.
    config
        FDR knobs.

    Returns
    -------
    dict[str, DEResult]
        Same dict with the three new columns
        (``padj_within_celltype``, ``padj_global``, ``strong_evidence``)
        on each table. The function modifies tables in-place and returns
        the same dict for chaining.
    """
    import pandas as pd

    cfg = config if config is not None else FDRConfig()

    # Pool every (cell_type, gene, pvalue) into one global vector
    pooled_rows: list[tuple[str, str, float]] = []
    for ct, res in de_results.items():
        for gene, p in zip(res.table.index, res.table[cfg.pvalue_col], strict=False):
            pooled_rows.append((ct, str(gene), float(p)))

    pooled = pd.DataFrame(pooled_rows, columns=["cell_type", "gene", "pvalue"])
    pooled["padj_global"] = bh_fdr(pooled["pvalue"].to_numpy())

    # Within-cell-type FDR (also recompute, since pyDESeq2's "padj" may
    # already be filled but we want to be explicit about which alpha we used)
    for ct, res in de_results.items():
        table = res.table.copy()
        table[cfg.within_padj_col] = bh_fdr(table[cfg.pvalue_col].to_numpy())
        # Merge in global padj
        global_for_ct = pooled.loc[pooled["cell_type"] == ct, ["gene", "padj_global"]].set_index(
            "gene"
        )
        table = table.join(global_for_ct, how="left")
        # Strong-evidence flag: both pass alpha
        table[cfg.strong_evidence_col] = (
            (table[cfg.within_padj_col] < cfg.alpha) & (table["padj_global"] < cfg.alpha)
        ).fillna(False)
        res.table = table
        n_within = int((table[cfg.within_padj_col] < cfg.alpha).sum())
        n_global = int((table["padj_global"] < cfg.alpha).sum())
        n_strong = int(table[cfg.strong_evidence_col].sum())
        logger.info(
            "FDR %s: %d significant within-cell-type, %d global, %d strong.",
            ct,
            n_within,
            n_global,
            n_strong,
        )

    return de_results
