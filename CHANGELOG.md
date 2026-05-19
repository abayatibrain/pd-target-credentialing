# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — biology layer (2026-05-19)
- `tests/fixtures/toy_smajic.h5ad` — synthetic 500-nucleus AnnData
  (2 donors × 3 cell types × 613 genes incl. the ADR-0005 marker panel,
  PD-relevant genes, and housekeeping anchors). Deterministic generator
  at `tests/fixtures/make_toy_anndata.py`.
- `src/pd_target_credentialing/annotate/markers.py` — ADR-0005's
  canonical nigra marker panel as a typed, immutable Python structure
  (9 cell types). HGNC-validation helper.
- `src/pd_target_credentialing/qc/nuclei_qc.py` — ADR-0002 thresholds
  (500 genes / 5% mito / Scrublet). Per-sample retention reporting with
  <50% flagging per brief §4.7. Scrublet wrapped behind ImportError
  fallback so minimal CI passes without the heavy dep.
- `src/pd_target_credentialing/io/smajic2022.py` — Smajić 2022 loader.
  Toy mode reads the committed fixture; real mode (GEO download with
  SHA256 verification) is a `NotImplementedError` placeholder for the
  v0.2.0 milestone.
- `src/pd_target_credentialing/io/kamath2022.py` — Kamath 2022
  cross-cohort loader. Carries the SOX6+/CALB1+ subtype passthrough
  per ADR-0005 Q5.3.
- `src/pd_target_credentialing/annotate/celltypes.py` — marker-based
  annotation with 0.15 ambiguity-margin filter (ADR-0005 Q5.2).
- `src/pd_target_credentialing/annotate/da_subtypes.py` — DA-subtype
  cross-check API contract (ADR-0005 Q5.3). v1.0.0 production semantics
  require the real Kamath embedding; toy-positional alignment ships now
  for plumbing tests.
- `src/pd_target_credentialing/de/pseudobulk.py` — donor × cell-type
  pseudobulk aggregation with 10-nucleus minimum (ADR-0006 Q6.2). Auto-
  detects `ancestry_pc*` covariates per Q6.3.
- `src/pd_target_credentialing/de/pydeseq2_runner.py` — pyDESeq2 wrapper
  for the headline DE pipeline. Heavy dep is optional at the library
  layer; integration tests use `pytest.importorskip("pydeseq2")`.
- `src/pd_target_credentialing/de/fdr.py` — per-cell-type and global
  BH-FDR per ADR-0007, plus a `strong_evidence` flag (passes both).
- `tests/unit/test_biology_layer.py` — 22 new tests for the biology
  modules. Brings total to 58 tests (+1 pyDESeq2 importorskip).

### Changed — pre-commit + tooling
- `pyproject.toml` — added `[[tool.mypy.overrides]]` block ignoring
  missing stubs for scientific-Python libs (anndata, scanpy, scipy,
  scrublet, harmonypy, pydeseq2, pandas).
- `.gitignore` — exempted `tests/fixtures/**/*.h5ad` and `*.parquet`
  from the default data-file exclusion (small synthetic fixtures are
  allowed).
- `.github/workflows/ci.yml` — extended mypy and pytest scopes to
  cover the new biology modules.

### Test surface as of this entry
- 58 passing tests, 1 skipped (pydeseq2 absent in minimal CI image).
- Coverage 86.69% across implemented modules; every module exceeds the
  §2.4 70% floor.

### Added
- Initial repository scaffolding per the Cowork brief §2.1.
- CI workflow (lint + type + test + coverage), docs workflow, release workflow.
- mkdocs-material site skeleton.
- **ADR sprint covering all nine §4.4 reasoning checkpoints:**
  - ADR-0001 — Choice of primary substantia nigra dataset (Smajić 2022
    primary; Kamath 2022 cross-validation; Wang 2024 deferred).
  - ADR-0002 — Per-nucleus QC thresholds (500 genes/nucleus, 5% max mito,
    Scrublet for doublets).
  - ADR-0003 — Normalization (log1p with library-size scaling to 10⁴);
    primary DE operates on raw counts so this choice is non-critical.
  - ADR-0004 — Batch integration (Harmony, integrating over donor not
    disease).
  - ADR-0005 — Cell-type annotation (marker-based primary with explicit
    nigra panel; Kamath cross-check for DA neurons).
  - ADR-0006 — Differential expression (pseudobulk + pyDESeq2 per cell
    type; per-cell Wilcoxon as direction-of-effect check only). Cites
    Squair 2021 on Type-I-error inflation in per-cell DE.
  - ADR-0007 — Multiple testing (per-cell-type BH-FDR at α=0.05 primary;
    global FDR reported alongside).
  - ADR-0008 — Composite score formula with weights `(0.5 genetic,
    0.3 literature, 0.2 animal)` and a bounded DE bonus.
  - ADR-0009 — Pre-registered calibration check against anchor genes;
    conditionally supersedes ADR-0008 on calibration failure.
- Decision-log index with dependency graph and reader-by-role guidance.
- QUESTIONS.md reorganized as a per-ADR review queue for Saturday-morning
  review.
- STATUS.md updated for the ADR-sprint phase.

### Notes
- ~~ADRs 0001–0009 are in **Proposed — awaiting Armin sign-off**
  status.~~ **Updated 2026-05-18: all nine ADRs ratified by Armin. See
  the sign-off entry below.** Biology-touching code implementation is
  now unblocked.

### ADR sign-off (Armin 2026-05-18)
All nine biology / methodology ADRs flipped from **Proposed** to
**Accepted**. Full response trail recorded inline in `QUESTIONS.md`
(27 timestamped responses across Q1.1–Q9.4 and three cross-cutting
items Q-X1/X2/X3). Substantive content changes from sign-off:

- **ADR-0005**: DA subtypes (SOX6+ vulnerable vs CALB1+ less-vulnerable)
  now **in scope for v1.0.0** (Q5.3).
- **ADR-0006**: Donor genetic-ancestry PCs added to the covariate set
  where dataset metadata supports them (Q6.3).
- **ADR-0008**: Confirmed the genetic/literature/animal **triad-only**
  scope for v1.0.0 — no OT drug or pathway channels added (Q8.3).
  Cross-disease panel expanded to include ≥1 oncology and ≥1
  inflammatory indication alongside AD/ALS/HD/FTD (Q-X1). All five
  tractability axes shown (Q-X2).
- **ADR-0009**: Positive anchor set expanded to **`{SNCA, GBA1, LRRK2,
  PRKN, PINK1, VPS35, PARK7}`** (Q9.1). PARK7 is the HGNC approved
  symbol for the historical DJ-1.
- **HGNC alias substitutions**: now also appear as a footnote in every
  rendered dossier (Q-X3), in addition to the WARNING log entry.

Cross-references to these clarifications are inlined in the relevant
ADR bodies under "Consequences" or alongside the original "Decision"
sections, with `(Armin sign-off 2026-05-18)` tags.

The decisions index (`docs/decisions/index.md`) and `STATUS.md` have
been updated to reflect the new phase.

### Added — foundational modules (engineering, §1.3 Cowork-decides)
- ADR-0010 — HGNC alias-resolver design (cache format, multi-mapping
  policy, alias-substitution logging).
- ADR-0011 — OpenTargets GraphQL client design (retry, cache,
  version pinning).
- ADR-0012 — Reactome ContentService client design (shares `_http/`
  helpers with OpenTargets).
- `src/pd_target_credentialing/_http/` — shared HTTP infrastructure:
  - `cache.py` — content-addressed on-disk cache keyed by
    `SHA256(method + url + canonicalized-body)`. Server-header-independent
    by design.
  - `retry.py` — bounded tenacity retry (3 attempts, 1s/2s/4s) with
    structured logging and a clean `HTTPRetryError` on exhaustion.
- `src/pd_target_credentialing/io/hgnc.py` — `HGNCResolver` with on-disk
  cache, `resolve` / `resolve_strict` / `resolve_many` APIs, alias
  substitution logged at WARNING, multi-mappings surfaced to a CSV for
  review (never silently merged).
- `src/pd_target_credentialing/evidence/opentargets.py` — `OpenTargetsClient`
  with pinned platform version, content-addressed cache, retry,
  and Pydantic response models. **Deliberately does not implement
  scoring** — that's ADR-0008 (pending sign-off).
- `src/pd_target_credentialing/evidence/reactome.py` — `ReactomeClient`
  returning HGNC-approved symbols (via the resolver) for pathway
  participants. Convenience constants for the two PD-relevant pathway
  IDs.
- `tests/unit/test_{http_helpers,hgnc,opentargets,reactome}.py` — 34
  tests, all HTTP traffic mocked via respx, no live network. Coverage
  93% across the new modules (every module exceeds the §2.4 70% floor).
- `tests/fixtures/http/{hgnc,opentargets,reactome}/` — checked-in JSON
  fixtures used by the respx mocks.

### Changed
- `tests/conftest.py` — autouse global-seed fixture is now tolerant of
  numpy/torch being absent (HTTP-only tests don't need them).
- `src/pd_target_credentialing/evidence/opentargets.py` — does not pass
  `base_url` to `httpx.Client` because httpx appends a trailing slash
  on empty-path requests, which silently diverges from the cache-key
  URL.
