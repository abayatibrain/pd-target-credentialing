# Status — week of 2026-05-18
Repo: pd-target-credentialing
Phase: **Implementation unblocked** — biology ADRs signed off; foundational
modules landed; biology code can begin next session.

## Completed this week
- Initial scaffolding committed per §2.1 of the Cowork brief.
- **Biology / methodology ADR sprint complete** — all nine §4.4
  reasoning checkpoints (ADR-0001 through ADR-0009) drafted, refined via
  the §1.4 two-pass critique, **and signed off by Armin on 2026-05-18.**
  All are now in **Accepted** status with response trail in
  [`QUESTIONS.md`](../QUESTIONS.md).
- **Engineering ADRs complete and accepted** (§1.3 Cowork-decides):
  ADR-0010 (HGNC), ADR-0011 (OpenTargets), ADR-0012 (Reactome).
- **Foundational modules implemented:**
  - `_http/` shared infrastructure (cache + retry).
  - `io/hgnc.py` — HGNC alias resolver with on-disk cache and §3.4
    alias-substitution logging.
  - `evidence/opentargets.py` — versioned, cached, retried GraphQL
    client.
  - `evidence/reactome.py` — ContentService client, HGNC-clean
    participant lists.
- **Tests:** 34 tests, 100% mocked HTTP (respx), zero live calls.
  Coverage: **93% across the four foundational modules** (every module
  exceeds §2.4's 70% floor).
- **CI green at commit `8258d90`** — lint, format, mypy, pytest with
  scoped 70% coverage floor all passing.

## ADRs Accepted this week (Armin 2026-05-18)
- ADR-0001 — Smajić 2022 primary, Kamath 2022 cross-cohort, Wang 2024
  deferred past v1.0.0.
- ADR-0002 — 500 genes/nucleus min, 5% max mito, Scrublet doublets.
- ADR-0003 — log1p with library-size scaling to 1e4.
- ADR-0004 — Harmony, integrating over donor (not disease).
- ADR-0005 — marker-based primary annotation with Kamath DA cross-check;
  DA subtypes (SOX6+/CALB1+) **in scope for v1.0.0**.
- ADR-0006 — pseudobulk + pyDESeq2; covariates = condition + age + sex
  + PMI + **donor genetic-ancestry PCs where metadata supports them**.
- ADR-0007 — per-cell-type BH-FDR α=0.05 primary, global FDR reported
  alongside.
- ADR-0008 — composite score `(0.5 · genetic + 0.3 · literature + 0.2 ·
  animal + bounded DE bonus ≤0.15)`. **Triad-only for v1.0.0** (no OT
  drug or pathway channels). Cross-disease panel includes
  non-neurodegenerative indications. All 5 tractability axes shown.
- ADR-0009 — calibration with **expanded positive anchor set** `{SNCA,
  GBA1, LRRK2, PRKN, PINK1, VPS35, PARK7}`; negative anchors
  `{ACTB, GAPDH, HPRT1, RPL13A, UBC}`. Cohen's d ≥ 1.0 + CI ≥ 0.5 +
  zero rank overlap is the pre-registered pass. Calibration failure
  halts the v1.0.0 tag.

## Cross-cutting clarifications (Armin 2026-05-18)
- Q-X1: cross-disease panel includes ≥1 oncology and ≥1 inflammatory
  indication alongside the neurodegenerative comparators (AD, ALS, HD,
  FTD).
- Q-X2: all five tractability axes shown — small molecule, antibody,
  PROTAC, ASO, gene therapy.
- Q-X3: alias substitution log appears as a footnote in every dossier.

## Blockers and questions for Armin
- **None blocking.** All decision-required items from QUESTIONS.md are
  resolved. The biology implementation work is unblocked.

## Plan for next week (now unblocked)
- Implement Smajić 2022 loader against a 500-nucleus toy fixture
  committed under `tests/fixtures/`.
- Implement Kamath 2022 cross-cohort loader (skeleton with subtype-label
  passthrough).
- Wire HGNC + OpenTargets + Reactome into a thin "evidence gather"
  module against SNCA as a sample target; render a minimal dossier
  card (no scoring yet).
- Begin the demo notebook (`notebooks/01_demo.ipynb`) wiring the
  end-to-end shape end-to-end on toy data.
- Implement the DA-subtype cross-check against Kamath (per Q5.3) as a
  v0.2.0 milestone.

## Burn rate
- Hours this week: ~12 (scaffolding + biology ADR sprint + engineering
  ADRs + foundational modules + tests + CI debugging + sign-off
  processing).
- Hours to `v0.1.0`: estimated 10-15 (Smajić + Kamath loaders, SNCA
  evidence gather, demo notebook stub on toy data).
