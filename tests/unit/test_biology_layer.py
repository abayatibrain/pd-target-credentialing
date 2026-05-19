"""End-to-end tests for the biology layer.

Exercises every new module against the toy fixture:

- markers: panel shape, HGNC contract
- nuclei_qc: stages drop the right things
- io.smajic2022 / io.kamath2022: toy-mode loads, schema verified
- annotate.celltypes: known cell types recovered, ambiguous flagged
- annotate.da_subtypes: cross-check returns labelled assignments
- de.pseudobulk: aggregation respects the 10-nucleus minimum, ancestry-PC
  fallback works
- de.fdr: BH math correct on synthetic p-values

The pyDESeq2 fit itself is exercised only when the library is importable;
otherwise its test is skipped with importorskip so the lean CI image
doesn't have to install scientific-stack deps.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from pd_target_credentialing.annotate.celltypes import AnnotationConfig, annotate_celltypes
from pd_target_credentialing.annotate.da_subtypes import assign_da_subtypes
from pd_target_credentialing.annotate.markers import (
    get_panel,
    panel_genes_present_in,
)
from pd_target_credentialing.de.fdr import apply_fdr, bh_fdr
from pd_target_credentialing.de.pseudobulk import (
    PseudobulkConfig,
    aggregate_pseudobulks,
)
from pd_target_credentialing.io.kamath2022 import load_kamath2022
from pd_target_credentialing.io.smajic2022 import (
    SmajicLoadError,
    load_smajic2022,
)
from pd_target_credentialing.qc.nuclei_qc import QCConfig, apply_qc


@pytest.fixture(scope="module")
def toy_adata() -> ad.AnnData:
    """Load the committed toy fixture once per test module."""
    return load_smajic2022(mode="toy")


# ----------------------------------------------------------------------
# Marker panel
# ----------------------------------------------------------------------


def test_marker_panel_has_expected_cell_types() -> None:
    panel = get_panel()
    expected = {
        "DA_neuron",
        "GABAergic_neuron",
        "glutamatergic_neuron",
        "astrocyte",
        "microglia",
        "oligodendrocyte",
        "OPC",
        "endothelial",
        "pericyte",
    }
    assert set(panel.cell_types) == expected


def test_marker_panel_is_immutable() -> None:
    panel = get_panel()
    entry = panel["DA_neuron"]
    assert entry.markers == ("TH", "SLC6A3", "DDC", "NR4A2", "KCNJ6")
    assert isinstance(panel.all_markers, tuple)


def test_panel_genes_present_in_toy_fixture(toy_adata: ad.AnnData) -> None:
    presence = panel_genes_present_in(list(toy_adata.var.index))
    # DA, astrocyte, oligo markers were intentionally included in the fixture
    da_markers = get_panel()["DA_neuron"].markers
    astro_markers = get_panel()["astrocyte"].markers
    oligo_markers = get_panel()["oligodendrocyte"].markers
    for m in (*da_markers, *astro_markers, *oligo_markers):
        assert presence[m] is True, f"expected {m} in toy fixture"


# ----------------------------------------------------------------------
# QC
# ----------------------------------------------------------------------


# For the toy fixture (~600 genes total, sparse Poisson background) the
# default 500-gene threshold from ADR-0002 drops everything. The threshold
# is correct for real Smajić data (15-20k expressed genes per nucleus);
# the tests below override it with a toy-scale value to test the *logic*
# at every stage.
_TOY_QC = QCConfig(min_genes_per_nucleus=50)


def test_qc_drops_high_mito_cells(toy_adata: ad.AnnData) -> None:
    filtered, report = apply_qc(toy_adata, _TOY_QC)
    assert report.n_input == toy_adata.n_obs
    assert filtered.n_obs <= toy_adata.n_obs
    # The fixture programs 20 high-mito "bad" nuclei; QC should drop at
    # least some of them under the 5% threshold
    assert report.n_dropped_max_mito > 0


def test_qc_reports_per_sample_retention(toy_adata: ad.AnnData) -> None:
    _, report = apply_qc(toy_adata, _TOY_QC)
    assert "D01" in report.per_sample_retention
    assert "D02" in report.per_sample_retention
    assert all(0.0 <= v <= 1.0 for v in report.per_sample_retention.values())


def test_qc_respects_custom_min_genes(toy_adata: ad.AnnData) -> None:
    # Crank the threshold up — should drop many more cells
    strict = apply_qc(toy_adata, QCConfig(min_genes_per_nucleus=10_000))
    lenient = apply_qc(toy_adata, QCConfig(min_genes_per_nucleus=1))
    assert strict[0].n_obs < lenient[0].n_obs


# ----------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------


def test_smajic_toy_mode_loads(toy_adata: ad.AnnData) -> None:
    assert toy_adata.n_obs == 500
    for col in ("donor_id", "condition", "age", "sex", "PMI"):
        assert col in toy_adata.obs.columns
    assert "is_mito" in toy_adata.var.columns


def test_smajic_real_mode_not_implemented() -> None:
    with pytest.raises(NotImplementedError):
        load_smajic2022(mode="real")


def test_smajic_unknown_mode_raises() -> None:
    with pytest.raises(ValueError):
        load_smajic2022(mode="cosmic")  # type: ignore[arg-type]


def test_smajic_missing_fixture_raises(tmp_path: Path) -> None:
    bogus = tmp_path / "no_such.h5ad"
    with pytest.raises(SmajicLoadError):
        load_smajic2022(mode="toy", fixture_path=bogus)


def test_kamath_toy_mode_returns_subtype_column() -> None:
    adata = load_kamath2022(mode="toy")
    assert "da_subtype" in adata.obs.columns
    seen = set(adata.obs["da_subtype"].unique())
    # toy mode assigns DA cells to SOX6_pos / CALB1_pos; non-DA to unresolved
    assert seen & {"SOX6_pos", "CALB1_pos", "unresolved"}


# ----------------------------------------------------------------------
# Annotation
# ----------------------------------------------------------------------


def test_annotation_recovers_known_cell_types(toy_adata: ad.AnnData) -> None:
    filtered, _ = apply_qc(toy_adata)
    annotated = annotate_celltypes(filtered)
    assert "celltype" in annotated.obs.columns
    assert "celltype_confidence" in annotated.obs.columns
    assert "celltype_ambiguous" in annotated.obs.columns
    # Among non-ambiguous cells, the fixture's true DA neurons should be
    # called DA_neuron at high rate (the fixture's marker programming is
    # strong by design)
    non_amb = annotated.obs[~annotated.obs["celltype_ambiguous"]]
    da_truth_mask = non_amb["celltype_true"] == "DA_neuron"
    da_called_mask = non_amb["celltype"] == "DA_neuron"
    if da_truth_mask.any():
        recall = (da_truth_mask & da_called_mask).sum() / da_truth_mask.sum()
        assert recall > 0.8, f"DA-neuron recall too low: {recall:.2f}"


def test_annotation_flags_ambiguous_cells(toy_adata: ad.AnnData) -> None:
    annotated = annotate_celltypes(toy_adata)
    # The fixture programs ~50 ambiguous cells with mixed DA+astro markers
    # — at least some should land in the ambiguous bucket
    n_amb = int(annotated.obs["celltype_ambiguous"].sum())
    assert n_amb > 0


def test_annotation_margin_threshold_changes_outcome(toy_adata: ad.AnnData) -> None:
    strict = annotate_celltypes(toy_adata, config=AnnotationConfig(ambiguity_margin=0.5))
    loose = annotate_celltypes(toy_adata, config=AnnotationConfig(ambiguity_margin=0.0))
    # A higher margin should flag more cells as ambiguous
    assert strict.obs["celltype_ambiguous"].sum() >= loose.obs["celltype_ambiguous"].sum()


# ----------------------------------------------------------------------
# DA-subtype cross-check
# ----------------------------------------------------------------------


def test_da_subtype_cross_check_returns_only_da(toy_adata: ad.AnnData) -> None:
    annotated = annotate_celltypes(toy_adata)
    kamath = load_kamath2022(mode="toy")
    assignments = assign_da_subtypes(annotated, kamath)
    # Every assignment should correspond to a DA neuron on the Smajic side
    assert all(a.smajic_label == "DA_neuron" for a in assignments)
    assert len(assignments) > 0


# ----------------------------------------------------------------------
# Pseudobulk
# ----------------------------------------------------------------------


def test_pseudobulk_aggregates_per_donor_celltype(toy_adata: ad.AnnData) -> None:
    filtered, _ = apply_qc(toy_adata, _TOY_QC)
    annotated = annotate_celltypes(filtered)
    pseudobulks = aggregate_pseudobulks(annotated)
    assert "DA_neuron" in pseudobulks
    # Two donors in the fixture, both should contribute > 10 nuclei to DA neurons
    da = pseudobulks["DA_neuron"]
    assert da.n_donors == 2
    assert all(n >= 10 for n in da.n_nuclei_per_donor.values())


def test_pseudobulk_excludes_ambiguous(toy_adata: ad.AnnData) -> None:
    annotated = annotate_celltypes(toy_adata)
    pseudobulks = aggregate_pseudobulks(annotated)
    # "ambiguous" must never be a key in the pseudobulks dict
    assert "ambiguous" not in pseudobulks


def test_pseudobulk_min_nuclei_filter(toy_adata: ad.AnnData) -> None:
    annotated = annotate_celltypes(toy_adata)
    pseudobulks = aggregate_pseudobulks(
        annotated, PseudobulkConfig(min_nuclei_per_pseudobulk=10_000)
    )
    # No donor x celltype pair has 10k nuclei in the toy fixture
    assert len(pseudobulks) == 0


def test_pseudobulk_picks_up_ancestry_pcs(toy_adata: ad.AnnData) -> None:
    annotated = annotate_celltypes(toy_adata)
    annotated.obs["ancestry_pc1"] = 0.5
    annotated.obs["ancestry_pc2"] = -0.3
    pseudobulks = aggregate_pseudobulks(annotated)
    da_meta = pseudobulks["DA_neuron"].metadata
    assert "ancestry_pc1" in da_meta.columns
    assert "ancestry_pc2" in da_meta.columns


# ----------------------------------------------------------------------
# FDR
# ----------------------------------------------------------------------


def test_bh_fdr_basic() -> None:
    pvalues = np.array([0.01, 0.02, 0.5, 0.9])
    out = bh_fdr(pvalues)
    assert out.shape == pvalues.shape
    assert np.all(out >= pvalues - 1e-10)
    assert np.all(out <= 1.0)


def test_bh_fdr_monotonic_under_sort() -> None:
    pvalues = np.array([0.001, 0.01, 0.04, 0.05, 0.2, 0.8])
    out = bh_fdr(pvalues)
    sorted_out = out[np.argsort(pvalues)]
    assert np.all(np.diff(sorted_out) >= -1e-10)


def test_bh_fdr_handles_nan() -> None:
    pvalues = np.array([0.01, np.nan, 0.5])
    out = bh_fdr(pvalues)
    assert np.isnan(out[1])
    assert not np.isnan(out[0])


def test_bh_fdr_empty_input() -> None:
    out = bh_fdr(np.array([]))
    assert out.shape == (0,)


def test_apply_fdr_adds_both_columns() -> None:
    from pd_target_credentialing.de.pydeseq2_runner import DEResult

    table = pd.DataFrame(
        {"pvalue": [0.001, 0.01, 0.05, 0.5]},
        index=["g1", "g2", "g3", "g4"],
    )
    de = {"DA_neuron": DEResult(cell_type="DA_neuron", table=table, n_donors=2, n_genes_tested=4)}
    out = apply_fdr(de)
    table_out = out["DA_neuron"].table
    assert "padj_within_celltype" in table_out.columns
    assert "padj_global" in table_out.columns
    assert "strong_evidence" in table_out.columns


# ----------------------------------------------------------------------
# pyDESeq2 integration (skipped if heavy dep missing)
# ----------------------------------------------------------------------


def test_pydeseq2_smoke(toy_adata: ad.AnnData) -> None:
    pytest.importorskip("pydeseq2")
    from pd_target_credentialing.de.pydeseq2_runner import run_pydeseq2_per_celltype

    filtered, _ = apply_qc(toy_adata)
    annotated = annotate_celltypes(filtered)
    pseudobulks = aggregate_pseudobulks(annotated)
    results = run_pydeseq2_per_celltype(pseudobulks)
    assert "DA_neuron" in results
    da = results["DA_neuron"].table
    assert "log2FoldChange" in da.columns
    assert "pvalue" in da.columns
