# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
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
- ADRs 0001–0009 are in **Proposed — awaiting Armin sign-off** status.
  Biology-touching code implementation does not begin until sign-off
  lands in `QUESTIONS.md`.

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
