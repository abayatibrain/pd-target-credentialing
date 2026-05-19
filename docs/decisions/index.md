# Decision log

Every non-trivial decision in this repo lives here as an Architecture
Decision Record (ADR). Reading these in order should let you reconstruct
every meaningful judgment call that shaped the code, without reading the
code itself.

The ADR template lives at [`templates/adr.md`](../templates/adr.md).

## Index — current ADRs

### Biology / methodology (Armin-ratified 2026-05-18)

| ADR | Title | Status |
|-----|---------------------------------------------------------|--------|
| [0001](0001.md) | Choice of primary substantia nigra dataset | Accepted |
| [0002](0002.md) | Per-nucleus QC thresholds | Accepted |
| [0003](0003.md) | Normalization (log1p vs scran vs SCT) | Accepted |
| [0004](0004.md) | Batch integration (Harmony vs BBKNN vs scVI) | Accepted |
| [0005](0005.md) | Cell-type annotation strategy | Accepted |
| [0006](0006.md) | Differential expression method | Accepted |
| [0007](0007.md) | Multiple-testing correction | Accepted |
| [0008](0008.md) | OpenTargets evidence weighting and score formula | Accepted |
| [0009](0009.md) | Confidence-score calibration via anchor genes | Accepted |

### Engineering (Cowork-decides per §1.3)

| ADR | Title | Status |
|-----|---------------------------------------------------------|--------|
| [0010](0010.md) | HGNC alias-resolver design | Accepted |
| [0011](0011.md) | OpenTargets GraphQL client design | Accepted |
| [0012](0012.md) | Reactome ContentService client design | Accepted |

All twelve ADRs are now Accepted. The biology implementation phase begins
next session. The full Armin response trail lives in
[`../QUESTIONS.md`](../../QUESTIONS.md).

## Dependency graph

```
ADR-0001 (dataset)
   ├── ADR-0002 (QC; thresholds chosen for nuclei)
   ├── ADR-0004 (integration; scale-matched to ~11 donors)
   └── ADR-0005 (annotation; cross-validated against Kamath)
        └── ADR-0006 (DE; consumes cell-type labels)
             ├── ADR-0007 (FDR; corrects DE p-values)
             └── ADR-0008 (score; consumes log2FC + padj)
                  └── ADR-0009 (calibration; ratifies or invalidates 0008)

ADR-0003 (normalization) is referenced by ADR-0004/0005 but NOT by ADR-0006,
because pseudobulk DE operates on raw counts and is robust to within-cell
normalization choice.
```

Engineering-side: ADR-0010 (HGNC) is the bedrock of ADR-0011 (OpenTargets)
and ADR-0012 (Reactome) — both clients return HGNC-approved symbols by
routing through the resolver.

## How to read these ADRs

If you are a hiring manager scanning this repo for ten minutes, start with
ADR-0006 (DE method) and ADR-0009 (calibration). Those two ADRs encode the
most consequential biology and statistics decisions in the pipeline.

If you are a biology PI evaluating credibility, also read ADR-0001 and
ADR-0008.

If you are an engineer evaluating reproducibility, also read ADR-0003,
ADR-0004, and ADR-0007.
