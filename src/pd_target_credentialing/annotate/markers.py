"""Canonical marker gene table for substantia nigra cell types.

The panel below is the load-bearing artifact ratified by Armin in ADR-0005
(sign-off 2026-05-18, QUESTIONS.md Q5.1). Every cell-type annotation in
this pipeline scores nuclei against these markers; the panel itself is
the citation of record.

Symbols use HGNC-approved nomenclature. Where a historical/alias name is
likely to surface in raw data (e.g., ``DAT`` for ``SLC6A3``), it is *not*
listed here — the alias resolution happens through the
:class:`pd_target_credentialing.io.hgnc.HGNCResolver` at data ingest, not
inside the panel.

Citations
---------
- Kamath T. *et al.* (2022) *Nat Neurosci* 25:588-595 — DA subtype markers.
- Smajić S. *et al.* (2022) *Brain* 145(3):964-978 — primary nigra atlas
  whose annotations seeded the panel.
- Hu Y. *et al.* (2017) *Bioinformatics* 33(2):248-250 — marker
  databases consulted for non-neuronal types.

Example
-------
>>> panel = get_panel()
>>> "DA_neuron" in panel.cell_types
True
>>> panel["DA_neuron"].markers
('TH', 'SLC6A3', 'DDC', 'NR4A2', 'KCNJ6')
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator

    from pd_target_credentialing.io.hgnc import HGNCResolver


@dataclass(frozen=True)
class CellTypeMarkers:
    """Markers and provenance notes for a single nigra cell type."""

    cell_type: str
    """Internal label used throughout the pipeline (snake-case)."""

    markers: tuple[str, ...]
    """HGNC-approved gene symbols. Order is the panel's order; the
    annotation module does not assume any particular ordering."""

    notes: str
    """Free-form notes on biology / known caveats. Surfaces in dossier
    footnotes when the cell type is reported."""


@dataclass(frozen=True)
class MarkerPanel:
    """The full panel for nigra annotation. Read via :func:`get_panel`."""

    entries: tuple[CellTypeMarkers, ...]

    @property
    def cell_types(self) -> tuple[str, ...]:
        """Cell-type labels in panel order."""
        return tuple(e.cell_type for e in self.entries)

    @property
    def all_markers(self) -> tuple[str, ...]:
        """Every marker across every cell type, deduplicated, in panel order."""
        seen: set[str] = set()
        out: list[str] = []
        for e in self.entries:
            for m in e.markers:
                if m not in seen:
                    seen.add(m)
                    out.append(m)
        return tuple(out)

    def __getitem__(self, cell_type: str) -> CellTypeMarkers:
        for e in self.entries:
            if e.cell_type == cell_type:
                return e
        raise KeyError(f"{cell_type!r} not in marker panel")

    def __iter__(self) -> Iterator[CellTypeMarkers]:
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)


# ----------------------------------------------------------------------
# The canonical panel — ratified ADR-0005, Armin 2026-05-18
# ----------------------------------------------------------------------

_PANEL = MarkerPanel(
    entries=(
        CellTypeMarkers(
            cell_type="DA_neuron",
            markers=("TH", "SLC6A3", "DDC", "NR4A2", "KCNJ6"),
            notes="KCNJ6 (Girk2) enriches in the vulnerable A9 subtype.",
        ),
        CellTypeMarkers(
            cell_type="GABAergic_neuron",
            markers=("GAD1", "GAD2", "SLC32A1"),
            notes="SLC32A1 is VGAT.",
        ),
        CellTypeMarkers(
            cell_type="glutamatergic_neuron",
            markers=("SLC17A7", "SLC17A6"),
            notes="VGLUT1 / VGLUT2.",
        ),
        CellTypeMarkers(
            cell_type="astrocyte",
            markers=("AQP4", "GFAP", "SLC1A2"),
            notes="SLC1A2 is EAAT2.",
        ),
        CellTypeMarkers(
            cell_type="microglia",
            markers=("CSF1R", "P2RY12", "TMEM119", "CX3CR1"),
            notes="Activated microglia downregulate P2RY12 and TMEM119.",
        ),
        CellTypeMarkers(
            cell_type="oligodendrocyte",
            markers=("PLP1", "MOG", "MBP"),
            notes="Mature OL panel; OPC markers are separate.",
        ),
        CellTypeMarkers(
            cell_type="OPC",
            markers=("PDGFRA", "CSPG4"),
            notes="Oligodendrocyte progenitor cells.",
        ),
        CellTypeMarkers(
            cell_type="endothelial",
            markers=("CLDN5", "PECAM1"),
            notes="CLDN5 enriches in brain endothelium.",
        ),
        CellTypeMarkers(
            cell_type="pericyte",
            markers=("PDGFRB", "RGS5"),
            notes="",
        ),
    )
)


def get_panel() -> MarkerPanel:
    """Return the canonical nigra marker panel.

    The returned object is immutable; share freely.

    Returns
    -------
    MarkerPanel
        The Armin-ratified panel from ADR-0005.
    """
    return _PANEL


def validate_panel_against_hgnc(
    resolver: HGNCResolver, *, panel: MarkerPanel | None = None
) -> dict[str, str]:
    """Sanity-check that every panel symbol is HGNC-approved.

    Parameters
    ----------
    resolver
        A :class:`HGNCResolver` (live or test-mocked).
    panel
        Panel to validate. Defaults to the canonical panel.

    Returns
    -------
    dict[str, str]
        Mapping from any *non-approved* input → the resolved approved symbol.
        An empty dict means the panel is clean. Multi-mappings and not-found
        results raise; the panel is by definition load-bearing and must not
        ship with ambiguous symbols.

    Raises
    ------
    KeyError
        If any marker resolves to "not found" against HGNC.
    ValueError
        If any marker is a multi-mapping.
    """
    panel = panel if panel is not None else _PANEL
    substitutions: dict[str, str] = {}
    for entry in panel:
        for symbol in entry.markers:
            # resolve_strict raises on multi-mapping or not-found; that's
            # what we want for a load-bearing panel.
            approved = resolver.resolve_strict(symbol)
            if approved != symbol:
                substitutions[symbol] = approved
    return substitutions


def panel_genes_present_in(var_index: list[str] | tuple[str, ...]) -> dict[str, bool]:
    """Check which panel markers exist in a dataset's gene index.

    Parameters
    ----------
    var_index
        Gene symbols from an AnnData's ``.var.index`` or equivalent.

    Returns
    -------
    dict[str, bool]
        Mapping from every panel marker → True if present in ``var_index``.
        Use this at load time to warn about missing markers before they
        cause silent annotation failures downstream.
    """
    present = set(var_index)
    return {m: m in present for m in _PANEL.all_markers}
