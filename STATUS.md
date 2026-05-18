# Status — week of 2026-05-18
Repo: pd-target-credentialing
Phase: Foundations landed; biology ADRs still awaiting Armin sign-off

## Completed this week
- Initial scaffolding committed per §2.1 of the Cowork brief.
- **Biology / methodology ADR sprint complete** — all nine §4.4
  reasoning checkpoints drafted as proper ADRs in `docs/decisions/`
  (ADR-0001 through ADR-0009). All in **Proposed — awaiting Armin
  sign-off** status.
- **Engineering ADRs complete and accepted** (§1.3 Cowork-decides):
  ADR-0010 (HGNC), ADR-0011 (OpenTargets), ADR-0012 (Reactome).
- **Foundational modules implemented:**
  - `_http/` shared infrastructure (cache + retry).
  - `io/hgnc.py` — HGNC alias resolver with on-disk cache and §3.4
    alias-substitution logging.
  - `evidence/opentargets.py` — versioned, cached, retried GraphQL
    client. Returns raw evidence; scoring deferred to ADR-0008 sign-off.
  - `evidence/reactome.py` — ContentService client, HGNC-clean
    participant lists.
- **Tests:** 34 tests, 100% mocked HTTP (respx), zero live calls.
  Coverage: **93% across the four foundational modules** (every
  module exceeds §2.4's 70% floor).

## ADRs added or updated this week
- ADR-0001 through ADR-0009: biology / methodology (Proposed).
- ADR-0010 through ADR-0012: engineering (Accepted).

## Blockers and questions for Armin
- `QUESTIONS.md` Q1.x through Q9.x — the Saturday-morning review
  queue. None of the foundational modules above touch any pending
  decision, so they can land independently. Once Q5.1 (marker panel),
  Q6.1–6.3 (DE method), and Q8/Q9 (score formula + calibration) are
  ratified, biology-touching code is unblocked.

## Plan for next week (assuming Armin replies land)
- Implement Smajić 2022 loader against a 500-nucleus toy fixture
  committed under `tests/fixtures/`.
- Wire HGNC + OpenTargets + Reactome into a thin "evidence gather"
  module against a sample target (SNCA) — no scoring yet, just the
  data-acquisition shape.
- Begin the demo notebook (`notebooks/01_demo.ipynb`).

## Burn rate
- Hours this week: ~10 (scaffolding + biology ADR sprint + engineering
  ADRs + foundational modules + tests).
- Hours to `v0.1.0`: estimated 8-12 once biology ADRs are signed off
  (Smajić loader + integration + the first end-to-end demo notebook
  on toy fixtures).
