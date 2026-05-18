# Biology primer for pd-target-credentialing

Audience: ML / engineering readers who need the biological context to read
this repo's README confidently. Skip if you already know the territory.

## Parkinson's disease in one paragraph
Parkinson's disease is a synucleinopathy: the protein α-synuclein (gene
**SNCA**) aggregates into intracellular Lewy bodies, and the dopaminergic
neurons of the substantia nigra pars compacta die. The clinical syndrome
that follows is bradykinesia, rigidity, tremor, and postural instability,
accompanied by autonomic, cognitive, sleep, and olfactory features.

## The genetics map (minimum useful version)
| Gene | Inheritance | Mechanism (very abbreviated) |
|---|---|---|
| SNCA | AD; GWAS | α-synuclein aggregation |
| LRRK2 | AD (G2019S common) | Kinase gain-of-function; lysosomal |
| GBA1 | Heterozygous risk | Glucocerebrosidase deficiency; lysosomal |
| PRKN (Parkin) | AR | E3 ubiquitin ligase; mitophagy |
| PINK1 | AR | Mitochondrial kinase; mitophagy initiator |
| VPS35 | AD | Retromer; endosomal trafficking |
| DJ-1 (PARK7) | AR | Oxidative stress response |

The convergent biological themes are mitochondrial quality control,
lysosomal/autophagic function, and α-synuclein homeostasis. The lens this
repo uses is therapeutic target credentialing — for any candidate gene,
can we converge enough public evidence to defend a drug program?

## Why pseudobulk DE
Per-cell statistical tests treat each cell as an independent sample, which
overstates power by ignoring within-subject correlation. Aggregating to
subject × cell-type pseudobulks before running DESeq2 produces calibrated
subject-level inference. This repo uses pseudobulk DE for primary results
and reports per-cell Wilcoxon only as a secondary check.

## Authoritative sources
This primer was written from public, citable sources. Where a claim is made
about disease biology, the underlying source is one of:

- HGNC (https://www.genenames.org/) — gene symbols
- OpenTargets (https://platform.opentargets.org/) — target-disease associations
- Reactome (https://reactome.org/) — pathway definitions
- UniProt (https://www.uniprot.org/) — protein function
- Primary literature (cited in README §Method)

If you find a claim here that is not defensible from these sources, open
an issue — that is a defect.
