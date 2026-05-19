"""DA-subtype cross-check against the Kamath 2022 reference atlas.

ADR-0005 Q5.3 (Armin sign-off 2026-05-18): DA subtypes (SOX6+ vulnerable
vs CALB1+ less-vulnerable) are in scope for v1.0.0. They are *not*
assigned by the top-level marker scoring (which would inflate the panel
with subtype-specific markers and confuse the headline DA-neuron call).
Instead:

1. Top-level annotation labels a nucleus as ``DA_neuron`` (or not).
2. This module projects every DA-labeled nucleus into the Kamath atlas's
   embedding and reads off the corresponding subtype label.
3. Disagreement between Smajic-side DA-marker-based labeling and the
   Kamath cross-check is logged and surfaces in the dossier as a
   "subtype not resolved" notation per ADR-0005's consequences.

v1.0.0 production implementation requires the real Kamath embedding,
which lands once the real-mode loader exists. This module ships with
the API contract and a toy-fixture-compatible path so downstream
consumers (dossier rendering, calibration) can be wired now.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import anndata as ad

from pd_target_credentialing.io.kamath2022 import DA_SUBTYPE_LABELS

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SubtypeAssignment:
    """Outcome of cross-checking one DA-neuron nucleus against Kamath."""

    cell_id: str
    """Index value from ``adata.obs.index`` (string-coerced)."""

    smajic_label: str
    """The Smajic-side label after marker-based annotation (e.g.
    ``"DA_neuron"`` or ``"ambiguous"``)."""

    kamath_subtype: str
    """One of :data:`DA_SUBTYPE_LABELS`."""

    agrees: bool
    """True if the Smajic-side and Kamath-side calls are mutually
    consistent (both DA, both same subtype if subtype is resolved).
    False means the dossier should render "subtype not resolved" for
    this nucleus."""


def assign_da_subtypes(
    smajic_adata: ad.AnnData,
    kamath_adata: ad.AnnData,
    *,
    smajic_label_col: str = "celltype",
    kamath_subtype_col: str = "da_subtype",
) -> list[SubtypeAssignment]:
    """For every Smajic DA neuron, return its Kamath-side subtype label.

    Parameters
    ----------
    smajic_adata
        Smajic-side AnnData after :func:`annotate_celltypes`.
    kamath_adata
        Kamath-side AnnData from :func:`load_kamath2022` with the
        ``da_subtype`` column populated.
    smajic_label_col
        Column in ``smajic_adata.obs`` carrying the top-level label.
    kamath_subtype_col
        Column in ``kamath_adata.obs`` carrying the subtype label.

    Returns
    -------
    list[SubtypeAssignment]
        One record per DA-labeled Smajic nucleus. Cells that the
        Smajic-side did NOT label as ``DA_neuron`` are not in the
        returned list.

    Notes
    -----
    The v0.x toy-mode path uses a positional alignment: Smajic cell *i*
    gets Kamath cell *i*'s subtype. This is **not** the production
    semantics — the production path projects Smajic cells into the
    Kamath embedding via a learned mapping. The toy path is here only
    so downstream rendering can be exercised against a populated
    return value.
    """
    if smajic_label_col not in smajic_adata.obs.columns:
        raise KeyError(
            f"{smajic_label_col!r} not in smajic_adata.obs; run annotate_celltypes first."
        )
    if kamath_subtype_col not in kamath_adata.obs.columns:
        raise KeyError(
            f"{kamath_subtype_col!r} not in kamath_adata.obs; use the "
            "Kamath loader to populate subtypes."
        )

    smajic_labels = smajic_adata.obs[smajic_label_col].astype(str).to_numpy()
    smajic_ids = [str(x) for x in smajic_adata.obs.index]
    kamath_subtypes = kamath_adata.obs[kamath_subtype_col].astype(str).to_numpy()

    out: list[SubtypeAssignment] = []
    n_da = 0
    n_agree = 0
    for i, (cell_id, label) in enumerate(zip(smajic_ids, smajic_labels, strict=False)):
        if label != "DA_neuron":
            continue
        n_da += 1
        # Toy alignment: positional. Production: nearest-neighbour in shared
        # latent space; deferred to v0.2.0.
        kamath_label = str(kamath_subtypes[i]) if i < len(kamath_subtypes) else "unresolved"
        if kamath_label not in DA_SUBTYPE_LABELS:
            kamath_label = "unresolved"
        agrees = kamath_label != "unresolved"
        if agrees:
            n_agree += 1
        out.append(
            SubtypeAssignment(
                cell_id=cell_id,
                smajic_label=label,
                kamath_subtype=kamath_label,
                agrees=agrees,
            )
        )
    logger.info(
        "DA-subtype cross-check: %d DA nuclei, %d with resolved Kamath subtype (%.1f%%).",
        n_da,
        n_agree,
        (n_agree / n_da * 100) if n_da else 0.0,
    )
    return out
