"""Generate a synthetic toy AnnData fixture for biology-module tests.

Run once to produce ``tests/fixtures/toy_smajic.h5ad``. The fixture is small
(<200KB) and committed to the repo. Re-running with the same seed produces
byte-identical output.

Design choices (all engineering, no biology decisions):

- **Shape.** 500 nuclei x ~2000 genes. Enough cells per (donor x cell-type)
  to exceed the 10-nucleus pseudobulk minimum from ADR-0006.
- **Donors.** 2 (one PD, one control). This is intentionally tiny — the
  fixture exists to verify *code paths*, not statistical power.
- **Cell types.** 3 known labels: dopaminergic neuron, astrocyte,
  oligodendrocyte. Each defined by canonical markers from the
  ADR-0005 panel. A small "ambiguous" cluster is included so the
  annotation module's 0.15-margin filter has something to flag.
- **Counts.** Poisson with per-cell-type means; one PD-relevant gene
  (SNCA) has a programmed disease effect in DA neurons so DE tests have
  a positive signal to find.
- **Metadata.** donor_id, condition (PD/control), age, sex, PMI — the
  ADR-0006 covariate set. ancestry PCs are deliberately absent so the
  fallback path (Q6.3) gets test coverage on real data.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 0xC0FFEE
N_NUCLEI = 500
# Keeping the gene universe tight: enough for the markers, PD-relevant genes,
# housekeeping anchors, and a small number of filler genes so neighbour-graph
# / clustering smoke tests have realistic dimensionality. Reduced from a
# larger size to keep the fixture under the pre-commit 1024 KB gate.
N_GENES = 600

# Gene panels (must be HGNC-approved symbols).
DA_MARKERS = ["TH", "SLC6A3", "DDC", "NR4A2", "KCNJ6"]
ASTRO_MARKERS = ["AQP4", "GFAP", "SLC1A2"]
OLIGO_MARKERS = ["PLP1", "MOG", "MBP"]
GABA_MARKERS = ["GAD1", "GAD2", "SLC32A1"]  # used to construct an ambiguous cluster
MICROGLIA_MARKERS = ["CSF1R", "P2RY12", "TMEM119"]
OPC_MARKERS = ["PDGFRA", "CSPG4"]
ENDO_MARKERS = ["CLDN5", "PECAM1"]
PERI_MARKERS = ["PDGFRB", "RGS5"]
GLUT_MARKERS = ["SLC17A7", "SLC17A6"]

# PD-credentialing-relevant genes that aren't markers but must appear so the
# downstream "evidence gather" smoke test has something to look up.
PD_GENES_OF_INTEREST = [
    "SNCA",
    "GBA1",
    "LRRK2",
    "PRKN",
    "PINK1",
    "VPS35",
    "PARK7",  # DJ-1
]

# Housekeeping anchors for the ADR-0009 calibration check.
HOUSEKEEPING = ["ACTB", "GAPDH", "HPRT1", "RPL13A", "UBC"]


def _build_gene_universe(rng: np.random.Generator) -> list[str]:
    """Combine all marker, PD, and housekeeping genes with synthetic filler genes
    until we reach N_GENES total. Filler genes get HGNC-like symbols
    (``GENE_0001`` ...). The order is deterministic for a given seed.
    """
    fixed = (
        DA_MARKERS
        + ASTRO_MARKERS
        + OLIGO_MARKERS
        + GABA_MARKERS
        + MICROGLIA_MARKERS
        + OPC_MARKERS
        + ENDO_MARKERS
        + PERI_MARKERS
        + GLUT_MARKERS
        + PD_GENES_OF_INTEREST
        + HOUSEKEEPING
    )
    # Deduplicate (some markers might overlap) while preserving order.
    seen: set[str] = set()
    ordered: list[str] = []
    for g in fixed:
        if g not in seen:
            seen.add(g)
            ordered.append(g)

    n_filler = N_GENES - len(ordered)
    filler = [f"FILLER_{i:04d}" for i in range(n_filler)]
    return ordered + filler


def _assign_cell_types(rng: np.random.Generator) -> np.ndarray:
    """Assign cell-type labels to 500 nuclei.

    Distribution chosen so that every (donor x cell-type) pseudobulk has
    >10 contributing nuclei (the ADR-0006 minimum), and every cell type
    has enough cells per donor for the QC module's per-sample retention
    check to make sense.
    """
    # 200 DA, 130 astrocyte, 120 oligo, 50 ambiguous (used to test the
    # 0.15-margin filter from ADR-0005)
    labels = (
        ["DA_neuron"] * 200 + ["astrocyte"] * 130 + ["oligodendrocyte"] * 120 + ["ambiguous"] * 50
    )
    arr = np.array(labels, dtype=object)
    rng.shuffle(arr)
    return arr


def _generate_counts(
    cell_types: np.ndarray,
    conditions: np.ndarray,
    gene_names: list[str],
    rng: np.random.Generator,
) -> np.ndarray:
    """Generate per-cell counts. Marker genes get high expression in their
    own cell type; the PD effect on SNCA is programmed in DA neurons of
    PD donors only."""
    n_cells = len(cell_types)
    n_genes = len(gene_names)
    counts = np.zeros((n_cells, n_genes), dtype=np.int32)

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Background baseline for every gene
    background = rng.poisson(0.5, size=(n_cells, n_genes))
    counts += background.astype(np.int32)

    # Marker-driven expression
    marker_map = {
        "DA_neuron": DA_MARKERS,
        "astrocyte": ASTRO_MARKERS,
        "oligodendrocyte": OLIGO_MARKERS,
        "ambiguous": DA_MARKERS[:2] + ASTRO_MARKERS[:1],  # mix → ambiguous
    }
    for ct, markers in marker_map.items():
        mask = cell_types == ct
        for m in markers:
            j = gene_to_idx.get(m)
            if j is None:
                continue
            # log-mean ~ 4-5; Poisson draws produce counts ~80-150 with variance
            mu = 100.0 if ct != "ambiguous" else 35.0
            counts[mask, j] += rng.poisson(mu, size=mask.sum()).astype(np.int32)

    # PD effect on SNCA in DA neurons of PD donors
    snca_idx = gene_to_idx["SNCA"]
    pd_da_mask = (cell_types == "DA_neuron") & (conditions == "PD")
    counts[pd_da_mask, snca_idx] += rng.poisson(50.0, size=pd_da_mask.sum()).astype(np.int32)
    ctrl_da_mask = (cell_types == "DA_neuron") & (conditions == "control")
    counts[ctrl_da_mask, snca_idx] += rng.poisson(20.0, size=ctrl_da_mask.sum()).astype(np.int32)

    # Housekeeping baseline elevated for all cells
    for h in HOUSEKEEPING:
        j = gene_to_idx[h]
        counts[:, j] += rng.poisson(30.0, size=n_cells).astype(np.int32)

    return counts


def make_toy() -> ad.AnnData:
    """Construct the toy AnnData."""
    rng = np.random.default_rng(SEED)

    gene_names = _build_gene_universe(rng)
    cell_types_true = _assign_cell_types(rng)

    # Two donors, balanced: 250 cells from each
    donor_ids = np.array(["D01"] * 250 + ["D02"] * 250, dtype=object)
    rng.shuffle(donor_ids)
    condition_map = {"D01": "PD", "D02": "control"}
    conditions = np.array([condition_map[d] for d in donor_ids], dtype=object)

    counts = _generate_counts(cell_types_true, conditions, gene_names, rng)

    # Add a mitochondrial-gene pattern so the 5% mito QC threshold has
    # something to filter against. Synthetic "MT-" prefixed genes appended.
    mito_genes = [f"MT-FILLER-{i:02d}" for i in range(13)]
    mito_counts = rng.poisson(2.0, size=(N_NUCLEI, len(mito_genes))).astype(np.int32)
    counts = np.hstack([counts, mito_counts])
    gene_names = gene_names + mito_genes
    # Add a few high-mito-content "bad nuclei" so the QC module has a real
    # signal to drop
    bad_idx = rng.choice(N_NUCLEI, size=20, replace=False)
    for j in range(len(mito_genes)):
        counts[bad_idx, len(gene_names) - len(mito_genes) + j] += rng.poisson(
            150.0, size=20
        ).astype(np.int32)

    # Build the AnnData
    obs = pd.DataFrame(
        {
            "donor_id": donor_ids,
            "condition": conditions,
            "age": np.where(donor_ids == "D01", 72.0, 68.0),
            "sex": np.where(donor_ids == "D01", "M", "F"),
            "PMI": np.where(donor_ids == "D01", 5.5, 7.2),
            "celltype_true": cell_types_true,  # ground truth label for tests
        },
        index=[f"cell_{i:04d}" for i in range(N_NUCLEI)],
    )
    var = pd.DataFrame(
        {
            "gene_name": gene_names,
            "is_mito": [g.startswith("MT-") for g in gene_names],
        },
        index=gene_names,
    )

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.uns["fixture_meta"] = {
        "seed": SEED,
        "schema_version": "1.0",
        "purpose": "synthetic-toy-for-pd-target-credentialing-tests",
        "do_not_redistribute_as_real_data": True,
    }
    return adata


def main() -> None:
    out = Path(__file__).parent / "toy_smajic.h5ad"
    adata = make_toy()
    adata.write_h5ad(out, compression="gzip")
    size_kb = out.stat().st_size // 1024
    print(f"wrote {out}  ({size_kb} KB, shape={adata.shape})")


if __name__ == "__main__":
    main()
