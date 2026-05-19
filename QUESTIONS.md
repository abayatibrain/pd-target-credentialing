# Open questions for Armin — pd-target-credentialing

This file is **append-only with respect to Armin's responses**. Cowork never
edits Armin's responses, ever. The framing of questions may be reorganized
as ADRs are written — the prior versions live in git history.

When Armin replies inline, prefix with a timestamp:
`> Armin (2026-05-24): ...`

The default protocol when a question is open: Cowork **does not** make the
decision unilaterally. See §1.3 of the brief for the authority matrix.

---

## How to use this file
The repository now has nine ADRs in `docs/decisions/`. Each ADR carries
its own "Open questions for Armin" section. This file is the **rolled-up
review queue** — every ADR question that needs sign-off appears here in
the same order as the ADRs. The intended workflow is one Saturday-morning
pass through this file, with replies inline.

ADRs with no remaining open questions after your replies will move from
"Proposed" to "Accepted" status.

---

## ADR-0001 — Primary substantia nigra dataset

**Q1.1** — Confirm Smajić 2022 as primary, or override.

> Armin (2026-05-18): Confirmed Smajić 2022 as primary.

**Q1.2** — Confirm Kamath 2022 as the cross-cohort validator for
DA-anchored claims.

> Armin (2026-05-18): Confirmed Kamath 2022 as cross-cohort validator.

**Q1.3** — Confirm Wang 2024 is acceptable to defer past `v1.0.0`.

> Armin (2026-05-18): Confirmed — Wang 2024 deferred past v1.0.0.

---

## ADR-0002 — Per-nucleus QC thresholds

**Q2.1** — Confirm 500 genes/nucleus minimum, 5% max mitochondrial
percentage, Scrublet for doublets — or override any of the three.

> Armin (2026-05-18): Confirmed: 500 genes/nucleus min, 5% max mito, Scrublet.

**Q2.2** — Should the QC report be checked into the repo as a static PDF
per dataset version, or generated on the fly by the demo notebook?

> Armin (2026-05-18): Generated on-the-fly by the demo notebook.

---

## ADR-0003 — Normalization

**Q3.1** — Confirm log1p with library-size scaling to 10⁴, or override to
scran or SCT (acknowledging the cost of cascading re-runs of ADR-0004
and ADR-0005).

> Armin (2026-05-18): Confirmed — log1p with library-size scaling to 1e4.

---

## ADR-0004 — Batch integration

**Q4.1** — Confirm Harmony, or override to scVI / BBKNN with the
GPU/non-determinism cost acknowledged.

> Armin (2026-05-18): Harmony confirmed.

**Q4.2** — Confirm integration is over **donor** and *not* over disease
status (the latter would erase the signal we are trying to measure).

> Armin (2026-05-18): Integrate over donor, not over disease status. Confirmed.

---

## ADR-0005 — Cell-type annotation

**Q5.1** — Confirm or revise the canonical nigra marker panel listed in
ADR-0005. This is the single most important biology call in this ADR.

> Armin (2026-05-18): Marker panel as drafted in ADR-0005 is fine — no additions.

**Q5.2** — Confirm the 0.15 ambiguity-margin threshold or specify another.

> Armin (2026-05-18): 0.15 ambiguity margin is fine.

**Q5.3** — Should DA subtypes (SOX6+ vulnerable vs CALB1+ less-vulnerable)
be in scope for `v1.0.0`, or deferred?

> Armin (2026-05-18): DA subtypes IN SCOPE for v1.0.0 (Kamath cross-check on SOX6+/CALB1+ subtypes).

---

## ADR-0006 — Differential expression method

**Q6.1** — Confirm pseudobulk + pyDESeq2 as primary, or override.

> Armin (2026-05-18): Pseudobulk + pyDESeq2 confirmed as primary.

**Q6.2** — Confirm the 10-nuclei minimum per (donor × cell type) for
pseudobulk inclusion, or specify another.

> Armin (2026-05-18): 10-nuclei minimum confirmed.

**Q6.3** — Confirm the covariate set: condition + age + sex + PMI. Should
donor genetic-ancestry PCs also enter the model where metadata supports
them?

> Armin (2026-05-18): Confirmed: condition + age + sex + PMI. Include donor genetic-ancestry PCs where metadata supports them.

---

## ADR-0007 — Multiple-testing correction

**Q7.1** — Confirm α = 0.05, or override (some teams use α = 0.10 for
target-ID purposes — defensible either way).

> Armin (2026-05-18): α = 0.05 confirmed.

**Q7.2** — Confirm reporting both per-cell-type and global FDR in the
dossier, or just per-cell-type.

> Armin (2026-05-18): Yes — report both per-cell-type and global FDR in the dossier.

---

## ADR-0008 — OpenTargets evidence weighting

**Q8.1** — Confirm the weight triple `(0.5 genetic, 0.3 literature, 0.2
animal)`, or override.

> Armin (2026-05-18): Confirmed: (0.5 genetic, 0.3 literature, 0.2 animal).

**Q8.2** — Confirm the DE-bonus cap of 0.15 and the |log2FC| = 1.5
saturation threshold, or override.

> Armin (2026-05-18): Confirmed: DE bonus cap 0.15, |log2FC| = 1.5 saturation.

**Q8.3** — Should OT drug-evidence and pathway-evidence channels be
added to the formula, or is the genetic/literature/animal triad
sufficient for `v1.0.0`?

> Armin (2026-05-18): No — the current genetic/literature/animal triad is sufficient for v1.0.0.

**Q8.4** — Confirm the dossier framing: "score is an evidence-convergence
index, not a probability." Precise phrasing is your call.

> Armin (2026-05-18): Confirmed: framing is "score is an evidence-convergence index, not a probability."

---

## ADR-0009 — Calibration

**Q9.1** — Confirm the positive anchor set `{SNCA, GBA1, LRRK2, PRKN,
PINK1}`. Add VPS35? Add DJ-1 (PARK7)?

> Armin (2026-05-18): Add VPS35 and PARK7 (DJ-1) to positive anchors. Final set: {SNCA, GBA1, LRRK2, PRKN, PINK1, VPS35, PARK7}.

**Q9.2** — Confirm the negative anchor set `{ACTB, GAPDH, HPRT1, RPL13A,
UBC}`.

> Armin (2026-05-18): Confirmed: {ACTB, GAPDH, HPRT1, RPL13A, UBC}.

**Q9.3** — Confirm the pre-registered pass thresholds: Cohen's d ≥ 1.0,
CI lower bound ≥ 0.5, zero rank overlap between positive and negative
anchors. Override if too strict or too loose.

> Armin (2026-05-18): Confirmed: Cohen's d ≥ 1.0, CI lower bound ≥ 0.5, zero rank overlap.

**Q9.4** — Confirm that a failing calibration **halts the `v1.0.0` tag**
until ADR-0008 is revised and re-run.

> Armin (2026-05-18): Confirmed — a failing calibration halts the v1.0.0 tag.

---

## Cross-cutting (not tied to a single ADR)

**Q-X1** — **Cross-disease panel scope.** The dossier compares PD score
to scores in: AD, ALS, HD, FTD. Should we include any non-neurodegenerative
diseases (e.g., an oncology indication, IBD) to surface promiscuous
targets?

> Armin (2026-05-18): Include at least one non-neurodegenerative indication (e.g., oncology and an inflammatory disease) in the cross-disease panel to surface promiscuous targets.

**Q-X2** — **Tractability bucket display.** OpenTargets provides several
tractability axes (small molecule, antibody, PROTAC, ASO, gene therapy).
Show all five in the dossier, or filter to small molecule + ASO + gene
therapy as the most PD-relevant?

> Armin (2026-05-18): Show all five tractability axes.

**Q-X3** — **HGNC alias substitution logging.** When the pipeline ingests
data using a deprecated alias (e.g., PARK2 → PRKN), should the substitution
log appear in the dossier as a footnote, or only in the run log?

> Armin (2026-05-18): Yes — alias substitution log appears as a footnote in the dossier.
