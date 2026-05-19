# Status — week of 2026-05-19
Repo: pd-target-credentialing
Phase: **v0.2.0 biology layer landed** — annotation + DE pipeline ships
on a toy fixture. Next slice is the real-data download path + the demo
notebook + the score formula.

## Completed overnight (2026-05-19, autonomous)
- **Toy AnnData fixture** committed at `tests/fixtures/toy_smajic.h5ad`
  (~376 KB, 500 nuclei × 613 genes, 2 donors × 3 cell types). Generator
  script at `tests/fixtures/make_toy_anndata.py` so reviewers can
  reproduce the fixture deterministically.
- **`annotate/markers.py`** — Armin-ratified ADR-0005 panel as a typed,
  immutable Python data structure. 9 cell types, every symbol HGNC-clean.
- **`qc/nuclei_qc.py`** — ADR-0002 thresholds (500 genes / 5% mito /
  Scrublet). Per-sample retention rate reported; samples below 50%
  retention flagged for review rather than silently dropped (brief §4.7).
  Scrublet wrapped behind `try/except ImportError` so minimal CI images
  still pass.
- **`io/smajic2022.py` + `io/kamath2022.py`** — dual-mode loaders. Toy
  mode reads the committed fixture; real mode raises NotImplementedError
  with a clear pointer to the v0.2.0 download-implementation slice.
  Kamath loader carries the SOX6+/CALB1+ subtype passthrough per
  ADR-0005 Q5.3.
- **`annotate/celltypes.py`** — marker-score-per-cluster with the 0.15
  ambiguity-margin filter (ADR-0005 Q5.2). Ambiguous nuclei are flagged
  in `obs["celltype_ambiguous"]` and excluded from downstream DE per
  ADR-0005's consequences.
- **`annotate/da_subtypes.py`** — DA-subtype cross-check skeleton (ADR-0005
  Q5.3). v1.0.0 production semantics require the real Kamath embedding;
  this commit ships the API contract and a toy-positional-alignment
  path so downstream code is exercisable.
- **`de/pseudobulk.py`** — pseudobulk aggregation per (donor × cell
  type). 10-nuclei minimum per ADR-0006 Q6.2. Covariate set follows
  Q6.3: condition + age + sex + PMI + auto-detected ancestry PCs.
- **`de/pydeseq2_runner.py`** — pyDESeq2 wrapper. Heavy dep is optional
  at the library layer; tests that need a live fit use
  `pytest.importorskip("pydeseq2")`. The pipeline contract is fully
  testable without the dep installed.
- **`de/fdr.py`** — per-cell-type BH-FDR + global BH-FDR + a
  `strong_evidence` flag (passes both). Pure scipy/numpy.

## Test surface
- 58 tests pass, 1 skipped (pydeseq2 not in the minimal CI image).
- Coverage: **86.69%** across the implemented modules — every module
  exceeds the §2.4 70% floor (lowest is `qc/nuclei_qc.py` at 79%,
  highest is `_http/retry.py` at 100%).
- All HTTP traffic is mocked (respx). No live network in CI.

## ADRs added or updated
None this session — the 12 existing ADRs cover every decision made.

## Blockers and questions for Armin
- **None blocking.** Everything stayed within already-ratified ADR
  boundaries. No new biology decisions were made; engineering choices
  (e.g., toy-fixture design, pyDESeq2 importorskip pattern) are
  documented in module docstrings.

## Pending commit through GitHub Desktop
Approximately 17 modified / new files. Suggested commit message:
> `feat(biology): toy fixture + marker panel + QC + loaders + annotation + DE skeleton (ADRs 0002, 0005, 0006, 0007)`

Local CI gauntlet results pre-commit:
- ✓ ruff check
- ✓ ruff format --check
- ✓ mypy (15 source files)
- ✓ pytest (58 passed, 1 skipped, 86.69% coverage)

## Plan for the next session
- Real-mode download path for Smajić 2022 + Kamath 2022 (GEO/SCP
  fetcher with SHA256 verification, cache dir under
  `$XDG_CACHE_HOME/pd_target_credentialing/`).
- Wire the full pipeline end-to-end on the toy fixture in
  `notebooks/01_demo.ipynb` (smoke run; one cell per pipeline stage).
- Score formula (ADR-0008) — needs at least one real OT call to
  sanity-check the anchor-gene rankings before locking in.
- Calibration (ADR-0009) — depends on the score formula.
- Dossier HTML rendering — once the first SNCA pass produces real
  numbers.

## Burn rate
- Hours this session: ~3 (toy fixture + 7 modules + tests + verification)
- Hours to `v0.1.0`: ~4 (real-mode loaders + demo notebook + score
  formula on toy data)
