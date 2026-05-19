"""pyDESeq2 wrapper for per-cell-type DE.

ADR-0006 commits the headline DE to pseudobulk + pyDESeq2. This module
takes a :class:`PseudobulkResult` and runs pyDESeq2 against the design
formula ``~ condition + age + sex + PMI + <ancestry_pcs>``.

The pyDESeq2 dependency is heavy. If it is not available in the
environment, :func:`run_pydeseq2_per_celltype` raises ``ImportError``
with a clear message — callers can ``pytest.importorskip("pydeseq2")``
to skip integration tests in minimal CI images.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import pandas as pd

    from pd_target_credentialing.de.pseudobulk import PseudobulkResult

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class DEConfig:
    """Per-celltype DE run config (ADR-0006 / ADR-0007 defaults)."""

    design_factor: str = "condition"
    """The design variable (PD vs control)."""

    contrast: tuple[str, str, str] = ("condition", "PD", "control")
    """The shrinkage contrast passed to DeseqStats."""

    covariates: tuple[str, ...] = ("age", "sex", "PMI")
    """Additional covariates always included if present in metadata."""

    include_ancestry_pcs: bool = True
    """Pick up ``ancestry_pc*`` columns automatically if present
    (Q6.3 fallback semantics)."""

    refit_cooks: bool = True
    """pyDESeq2 outlier handling; default-on matches DESeq2 R behaviour."""


@dataclass
class DEResult:
    """One per cell type."""

    cell_type: str
    table: pd.DataFrame
    """Columns: gene, log2FoldChange, lfcSE, stat, pvalue, padj_within_celltype.
    Index: gene name. The global-FDR column is added by the FDR module.
    """
    n_donors: int
    n_genes_tested: int


def _require_pydeseq2() -> None:
    try:
        import pydeseq2  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "pydeseq2 is required for the headline DE pipeline. Install "
            "with `uv pip install pydeseq2`. Tests that don't require a "
            "live pyDESeq2 fit can `pytest.importorskip('pydeseq2')`."
        ) from exc


def run_pydeseq2_per_celltype(
    pseudobulks: dict[str, PseudobulkResult],
    *,
    config: DEConfig | None = None,
) -> dict[str, DEResult]:
    """Run pyDESeq2 per cell type and return per-gene DE tables.

    Parameters
    ----------
    pseudobulks
        Output of :func:`aggregate_pseudobulks`.
    config
        Run knobs.

    Returns
    -------
    dict[str, DEResult]
        Keyed by cell type. The per-cell-type FDR column
        ``padj_within_celltype`` is the within-cell-type BH-FDR per
        ADR-0007.

    Raises
    ------
    ImportError
        If pyDESeq2 isn't installed.
    """
    import re

    _require_pydeseq2()
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    cfg = config if config is not None else DEConfig()
    out: dict[str, DEResult] = {}

    for ct, pb in pseudobulks.items():
        meta = pb.metadata.copy()
        # Identify ancestry PC covariates that are actually present
        ancestry_cols = (
            [c for c in meta.columns if re.match(r"^ancestry_pc\d+$", c)]
            if cfg.include_ancestry_pcs
            else []
        )
        design_factors = [cfg.design_factor, *cfg.covariates, *ancestry_cols]
        design_factors = [f for f in design_factors if f in meta.columns]

        # pyDESeq2 expects counts as samples x genes
        counts = pb.counts.T.astype(int)
        meta = meta.reindex(counts.index)
        # Drop rows with NA in design factors
        meta_clean = meta.dropna(subset=design_factors)
        counts = counts.loc[meta_clean.index]

        if counts.shape[0] < 2:
            logger.warning(
                "Skipping DE for %s: only %d donors after NA filter.",
                ct,
                counts.shape[0],
            )
            continue

        dds = DeseqDataSet(
            counts=counts,
            metadata=meta_clean,
            design_factors=design_factors,
            refit_cooks=cfg.refit_cooks,
            quiet=True,
        )
        dds.deseq2()
        stats = DeseqStats(dds, contrast=list(cfg.contrast), quiet=True)
        stats.summary()
        results = stats.results_df.copy()
        # Standardise column names + rename padj for ADR-0007 semantics
        results = results.rename(columns={"padj": "padj_within_celltype"})
        out[ct] = DEResult(
            cell_type=ct,
            table=results,
            n_donors=int(counts.shape[0]),
            n_genes_tested=int(counts.shape[1]),
        )
        logger.info(
            "DE %s: %d donors x %d genes; %d genes with padj_within_celltype < 0.05.",
            ct,
            counts.shape[0],
            counts.shape[1],
            int((results["padj_within_celltype"] < 0.05).sum()),
        )
    return out
