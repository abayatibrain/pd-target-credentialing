"""Pseudobulk aggregation for cell-type-conditioned DE.

ADR-0006 (Armin sign-off 2026-05-18) commits to **pseudobulk + pyDESeq2**
as the primary DE method. This module handles the aggregation; the
pyDESeq2 fit lives in :mod:`pd_target_credentialing.de.pydeseq2_runner`.

Per Armin's responses (Q6.2, Q6.3):

- Pseudobulks with fewer than ``min_nuclei_per_pseudobulk`` (default 10)
  contributing nuclei are dropped.
- The sample-level metadata kept alongside the pseudobulks includes
  ``condition + age + sex + PMI`` plus any genetic-ancestry PCs that
  exist in ``adata.obs`` (columns matching ``ancestry_pc[0-9]+``).
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad
    import numpy as np
    import pandas as pd

logger = logging.getLogger(__name__)

ANCESTRY_PC_PATTERN = re.compile(r"^ancestry_pc\d+$")


@dataclass(frozen=True)
class PseudobulkConfig:
    """Knobs for the aggregation step (ADR-0006-ratified defaults)."""

    min_nuclei_per_pseudobulk: int = 10
    """Q6.2 — drop (donor, celltype) pairs below this."""

    donor_col: str = "donor_id"
    """Column carrying the per-donor grouping."""

    celltype_col: str = "celltype"
    """Column carrying the cell-type label (annotate_celltypes output)."""

    condition_col: str = "condition"
    """The design variable. PD vs control."""

    covariate_cols: tuple[str, ...] = ("age", "sex", "PMI")
    """Always-include covariates. Ancestry PCs are added dynamically if
    present in adata.obs (Q6.3)."""

    exclude_ambiguous: bool = True
    """When True, drop nuclei labeled "ambiguous" before aggregation
    (per ADR-0005's annotation contract)."""


@dataclass
class PseudobulkResult:
    """One per cell type."""

    cell_type: str
    counts: pd.DataFrame  # genes x donors
    metadata: pd.DataFrame  # donor-level covariates aligned to counts.columns
    n_nuclei_per_donor: dict[str, int]
    n_genes: int = field(init=False)
    n_donors: int = field(init=False)

    def __post_init__(self) -> None:
        # post_init can't assign on a frozen dataclass; this is not frozen.
        self.n_genes = self.counts.shape[0]
        self.n_donors = self.counts.shape[1]


def _ancestry_pc_columns(obs_columns: list[str]) -> list[str]:
    return sorted(c for c in obs_columns if ANCESTRY_PC_PATTERN.match(c))


def aggregate_pseudobulks(
    adata: ad.AnnData, config: PseudobulkConfig | None = None
) -> dict[str, PseudobulkResult]:
    """Aggregate raw counts to (donor x celltype) pseudobulks.

    Parameters
    ----------
    adata
        Post-annotation AnnData. Must have raw counts available at
        ``adata.layers["counts"]`` (the QC step preserves them).
    config
        Aggregation knobs.

    Returns
    -------
    dict[str, PseudobulkResult]
        Keyed by cell type. Each value carries the genes x donors count
        matrix and the donor-level metadata (aligned columns).

    Notes
    -----
    Donors with fewer than ``min_nuclei_per_pseudobulk`` contributing
    nuclei for a given cell type are dropped from that cell type's
    pseudobulk *only*; they may still appear in other cell types.
    """
    import pandas as pd
    import scipy.sparse as sp

    cfg = config if config is not None else PseudobulkConfig()

    if "counts" not in adata.layers:
        raise ValueError(
            "Pseudobulk aggregation requires raw counts at "
            "adata.layers['counts'] (the QC step preserves them)."
        )
    counts = adata.layers["counts"]
    if sp.issparse(counts):
        counts = counts.toarray()

    obs = adata.obs.copy()
    if cfg.exclude_ambiguous and "celltype_ambiguous" in obs.columns:
        keep_mask = ~obs["celltype_ambiguous"].astype(bool).to_numpy()
        obs = obs.loc[keep_mask]
        counts = counts[keep_mask, :]

    ancestry_cols = _ancestry_pc_columns(list(obs.columns))
    if ancestry_cols:
        logger.info(
            "Including %d ancestry PCs as covariates: %s",
            len(ancestry_cols),
            ancestry_cols,
        )
    else:
        logger.info(
            "No ancestry-PC columns in obs; falling back to "
            "condition + age + sex + PMI per ADR-0006 Q6.3."
        )
    metadata_cols = (
        cfg.condition_col,
        *cfg.covariate_cols,
        *ancestry_cols,
    )

    gene_names = list(adata.var.index)
    out: dict[str, PseudobulkResult] = {}

    for ct in sorted(obs[cfg.celltype_col].astype(str).unique()):
        if ct == "ambiguous":
            continue
        ct_mask = (obs[cfg.celltype_col].astype(str) == ct).to_numpy()
        ct_obs = obs.loc[ct_mask]
        ct_counts = counts[ct_mask, :]

        # Aggregate per donor
        donor_groups = ct_obs.groupby(cfg.donor_col, sort=True)
        donor_counts_cols: dict[str, np.ndarray] = {}
        donor_metadata_rows: list[dict[str, object]] = []
        n_per_donor: dict[str, int] = {}

        for donor, idx in donor_groups.indices.items():
            n_cells = len(idx)
            if n_cells < cfg.min_nuclei_per_pseudobulk:
                logger.info(
                    "Dropping pseudobulk (donor=%s, celltype=%s): only %d nuclei (< %d minimum).",
                    donor,
                    ct,
                    n_cells,
                    cfg.min_nuclei_per_pseudobulk,
                )
                continue
            donor_counts_cols[str(donor)] = ct_counts[idx, :].sum(axis=0)
            n_per_donor[str(donor)] = n_cells
            # Donor-level metadata: take the first row's values for each col
            row = ct_obs.iloc[idx[0]]
            donor_metadata_rows.append(
                {"donor_id": str(donor), **{c: row.get(c) for c in metadata_cols}}
            )

        if not donor_counts_cols:
            logger.warning(
                "Cell type %s has no donors meeting the min-nuclei threshold; skipping.",
                ct,
            )
            continue

        counts_df = pd.DataFrame(donor_counts_cols, index=gene_names)
        metadata_df = pd.DataFrame(donor_metadata_rows).set_index("donor_id")
        metadata_df = metadata_df.loc[counts_df.columns]

        out[ct] = PseudobulkResult(
            cell_type=ct,
            counts=counts_df,
            metadata=metadata_df,
            n_nuclei_per_donor=n_per_donor,
        )

    logger.info(
        "Pseudobulks built for %d cell types: %s",
        len(out),
        sorted(out.keys()),
    )
    return out
