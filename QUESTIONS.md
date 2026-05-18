# Open questions for Armin — pd-target-credentialing

This file is **append-only**. Cowork never edits Armin's responses. Use
timestamps. New questions go at the bottom under the next-empty heading.

The default protocol when a question is open: Cowork **does not** make the
decision unilaterally. See §1.3 of the brief for the authority matrix.

---

## Q1 — Dataset choice (ADR-0001) — *decision required*

ADR-0001 proposes Smajić 2022 as primary and Kamath 2022 as cross-cohort
validation, dropping Wang 2024 to a future extension. Confirm or override.

> Armin: <reply here, with timestamp>

## Q2 — Anchor genes for confidence-score calibration (§4.4 item 9) — *decision required*

Default proposal: GBA1 and LRRK2 as known positives (high human-genetics
evidence, tractable). ACTB as a known negative (housekeeping, no plausible
PD link). Are these the right anchors, or do you want a richer panel?

> Armin: <reply>

## Q3 — Confidence-score weighting (§4.4 item 8) — *decision required*

Proposal: weight OpenTargets evidence types as genetic > literature > animal
(e.g., 0.5 / 0.3 / 0.2). Pseudobulk DE log2FC contributes a sign-aware bonus
capped at +0.15. The final score is bootstrap-95% CI'd over evidence
subsampling. Defend this in ADR-0008 or override.

> Armin: <reply>

## Q4 — Cross-disease panel scope (§4.2 section 6) — *decision required*

The dossier compares PD score to scores in: AD, ALS, HD, FTD. Should we
include any non-neurodegenerative diseases to surface promiscuous targets?

> Armin: <reply>

## Q5 — Tractability bucket display — *Cowork can decide if Armin defers*

OpenTargets provides several tractability axes (small molecule, antibody,
PROTAC, ASO, gene therapy). Show all five, or filter to small molecule
+ ASO + gene therapy as the most PD-relevant?

> Armin: <reply (or "Cowork's call")>
