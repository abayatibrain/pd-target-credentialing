"""Per-nucleus QC filters for substantia nigra snRNA-seq.

Implements the thresholds ratified in ADR-0002 (Armin 2026-05-18):

- >= 500 genes detected per nucleus
- <= 5% mitochondrial transcript fraction (nuclei-appropriate; not whole-cell)
- Scrublet for doublet detection (per-sample)
- A sample retaining < 50% of its loaded nuclei is **flagged for manual
  review** rather than silently dropped.

The doublet step is optional at the library level: if Scrublet is not
installed (e.g., in a minimal CI image), the QC pass still runs with a
clear note in the report that doublet flagging was skipped. Production
runs should always have Scrublet available.

Example
-------
>>> from pd_target_credentialing.qc.nuclei_qc import apply_qc, QCConfig
>>> # filtered, report = apply_qc(adata, QCConfig())
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad
    import numpy as np

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class QCConfig:
    """Per-nucleus QC thresholds. Defaults are the ADR-0002 values.

    Departures from these defaults require a superseding ADR per the
    Cowork brief Section 1.2.
    """

    min_genes_per_nucleus: int = 500
    """ADR-0002 Q2.1. nuclei-appropriate (vs. whole-cell ~200)."""

    max_mito_fraction: float = 0.05
    """ADR-0002 Q2.1. nuclei-appropriate (vs. whole-cell ~0.20)."""

    expected_doublet_rate: float = 0.07
    """Used by Scrublet's per-sample threshold tuning. 7% is the
    canonical 10x Genomics-recommended value for loading densities
    typical of nigra snRNA-seq."""

    sample_retention_warn_threshold: float = 0.50
    """Below this retention rate the sample is flagged for manual review
    rather than silently dropped (the brief is explicit about this in
    Section 4.7 — silent over-filtering is a documented pitfall)."""

    mito_var_column: str = "is_mito"
    """Column in ``adata.var`` flagging mitochondrial genes. The Smajic
    and Kamath loaders both populate this; if it's missing we fall back
    to the ``MT-`` prefix heuristic on the gene index."""

    sample_obs_column: str = "donor_id"
    """Column in ``adata.obs`` identifying per-sample (per-donor)
    grouping. Used for per-sample retention reporting and per-sample
    Scrublet runs."""


@dataclass
class QCReport:
    """Audit-trail dict returned alongside the filtered AnnData."""

    n_input: int
    n_output: int
    n_dropped_min_genes: int
    n_dropped_max_mito: int
    n_dropped_doublet: int
    per_sample_retention: dict[str, float] = field(default_factory=dict)
    flagged_samples: list[str] = field(default_factory=list)
    doublet_method: str = "skipped"
    config: QCConfig | None = None

    @property
    def retention_rate(self) -> float:
        """Overall fraction of input nuclei retained."""
        if self.n_input == 0:
            return 0.0
        return self.n_output / self.n_input

    def as_dict(self) -> dict[str, object]:
        """Plain-dict view for JSON serialization / dossier rendering."""
        return {
            "n_input": self.n_input,
            "n_output": self.n_output,
            "retention_rate": self.retention_rate,
            "n_dropped_min_genes": self.n_dropped_min_genes,
            "n_dropped_max_mito": self.n_dropped_max_mito,
            "n_dropped_doublet": self.n_dropped_doublet,
            "per_sample_retention": dict(self.per_sample_retention),
            "flagged_samples": list(self.flagged_samples),
            "doublet_method": self.doublet_method,
        }


def _mito_mask(adata: ad.AnnData, config: QCConfig) -> np.ndarray:
    """Return a boolean array marking mitochondrial genes."""
    import numpy as np

    if config.mito_var_column in adata.var.columns:
        return adata.var[config.mito_var_column].to_numpy().astype(bool)
    # Fallback: pattern-match on the gene index.
    return np.array([str(g).startswith("MT-") for g in adata.var.index], dtype=bool)


def _per_nucleus_metrics(adata: ad.AnnData, mito_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (n_genes_per_nucleus, mito_fraction_per_nucleus)."""
    import numpy as np
    import scipy.sparse as sp

    X = adata.X
    if sp.issparse(X):
        n_genes = np.asarray((X > 0).sum(axis=1)).ravel()
        total_counts = np.asarray(X.sum(axis=1)).ravel()
        mito_counts = np.asarray(X[:, mito_mask].sum(axis=1)).ravel()
    else:
        Xa = np.asarray(X)
        n_genes = (Xa > 0).sum(axis=1)
        total_counts = Xa.sum(axis=1)
        mito_counts = Xa[:, mito_mask].sum(axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        mito_frac = np.where(total_counts > 0, mito_counts / total_counts, 0.0)
    return n_genes.astype(float), mito_frac.astype(float)


def _try_scrublet(adata: ad.AnnData, config: QCConfig) -> tuple[np.ndarray, str]:
    """Run Scrublet per sample. If the library isn't importable, return an
    all-False mask and doublet_method = "skipped"."""
    import numpy as np

    try:
        import scrublet as scr
    except ImportError:
        logger.warning(
            "Scrublet not installed; skipping doublet detection. "
            "The CI/production environment should include Scrublet "
            "(per ADR-0002 Q2.1)."
        )
        return np.zeros(adata.n_obs, dtype=bool), "skipped"

    is_doublet = np.zeros(adata.n_obs, dtype=bool)
    samples = adata.obs[config.sample_obs_column].astype(str)
    for sample in samples.unique():
        mask = (samples == sample).to_numpy()
        if mask.sum() < 30:
            # Scrublet needs a reasonable number of cells per sample.
            continue
        sub_X = adata[mask].X
        try:
            sim = scr.Scrublet(sub_X, expected_doublet_rate=config.expected_doublet_rate)
            _, calls = sim.scrub_doublets(verbose=False)
            is_doublet[mask] = np.asarray(calls, dtype=bool)
        except Exception as exc:  # pragma: no cover — defensive
            logger.warning(
                "Scrublet failed on sample %s: %s. Marking sample as non-doublet-filtered.",
                sample,
                exc,
            )
    return is_doublet, "scrublet"


def apply_qc(adata: ad.AnnData, config: QCConfig | None = None) -> tuple[ad.AnnData, QCReport]:
    """Filter nuclei and return the filtered AnnData + an audit report.

    Parameters
    ----------
    adata
        Input AnnData with raw counts in ``adata.X``.
    config
        Thresholds. Defaults to the ADR-0002 values.

    Returns
    -------
    tuple
        ``(filtered_adata, qc_report)``. The filtered AnnData is a copy
        (the input is never mutated). The report carries the audit trail
        the dossier needs to footnote what was filtered out.

    Notes
    -----
    Order of filters: min genes -> max mito -> doublet. Each is reported
    separately so the dossier can show what dropped out at each step.
    """
    cfg = config if config is not None else QCConfig()
    n_input = adata.n_obs
    mito_mask = _mito_mask(adata, cfg)
    n_genes, mito_frac = _per_nucleus_metrics(adata, mito_mask)

    # Stage 1: min genes
    pass_genes = n_genes >= cfg.min_genes_per_nucleus
    n_dropped_min_genes = int((~pass_genes).sum())

    # Stage 2: max mito (only evaluated on cells that passed stage 1)
    pass_mito = mito_frac <= cfg.max_mito_fraction
    n_dropped_max_mito = int((~pass_mito & pass_genes).sum())

    pre_doublet_mask = pass_genes & pass_mito
    pre_doublet = adata[pre_doublet_mask].copy()

    # Stage 3: doublets
    is_doublet, doublet_method = _try_scrublet(pre_doublet, cfg)
    n_dropped_doublet = int(is_doublet.sum())
    keep = ~is_doublet
    filtered = pre_doublet[keep].copy()

    # Per-sample retention
    per_sample_retention: dict[str, float] = {}
    flagged: list[str] = []
    if cfg.sample_obs_column in adata.obs.columns:
        for sample in adata.obs[cfg.sample_obs_column].astype(str).unique():
            in_count = int((adata.obs[cfg.sample_obs_column].astype(str) == sample).sum())
            out_count = int((filtered.obs[cfg.sample_obs_column].astype(str) == sample).sum())
            rate = (out_count / in_count) if in_count > 0 else 0.0
            per_sample_retention[sample] = rate
            if rate < cfg.sample_retention_warn_threshold:
                flagged.append(sample)
                logger.warning(
                    "Sample %s retained only %.1f%% of nuclei after QC — "
                    "flagged for manual review (threshold: %.1f%%).",
                    sample,
                    rate * 100,
                    cfg.sample_retention_warn_threshold * 100,
                )

    report = QCReport(
        n_input=n_input,
        n_output=filtered.n_obs,
        n_dropped_min_genes=n_dropped_min_genes,
        n_dropped_max_mito=n_dropped_max_mito,
        n_dropped_doublet=n_dropped_doublet,
        per_sample_retention=per_sample_retention,
        flagged_samples=flagged,
        doublet_method=doublet_method,
        config=cfg,
    )
    logger.info(
        "QC: %d -> %d nuclei retained (%.1f%%). "
        "Dropped: %d <min_genes, %d >mito, %d doublets (method=%s).",
        n_input,
        filtered.n_obs,
        report.retention_rate * 100,
        n_dropped_min_genes,
        n_dropped_max_mito,
        n_dropped_doublet,
        doublet_method,
    )
    return filtered, report
