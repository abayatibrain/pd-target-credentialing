# Data sources for pd-target-credentialing

No data files are committed. This document is the manifest;
`scripts/download_data.sh` is the executable form.

## Datasets

### Smajić 2022 — primary
- **Source**: GSE178265 (NCBI GEO) and ArrayExpress
- **Citation**: Smajić S. *et al.* (2022) *Brain* 145(3):964–978.
  doi:10.1093/brain/awab406
- **License**: CC-BY 4.0 on the count matrices via GEO
- **Approximate size**: ~2 GB raw
- **Snapshot date**: pinned at download time; recorded in
  `data/snapshots.json`

### Kamath 2022 — cross-cohort validation
- **Source**: Single Cell Portal SCP1768 (Broad Institute)
- **Citation**: Kamath T. *et al.* (2022) *Nat Neurosci* 25:588–595.
  doi:10.1038/s41593-022-01061-1
- **License**: dataset terms on the Single Cell Portal page
- **Approximate size**: ~6 GB raw

### OpenTargets — evidence layer
- **Source**: OpenTargets Platform GraphQL API
  (https://api.platform.opentargets.org/api/v4/graphql)
- **Citation**: Ochoa D. *et al.* (2023) *Nucleic Acids Research*
  51(D1):D1353-D1359. doi:10.1093/nar/gkac1046
- **Pinned release**: documented per-call; the local cache records the
  release version and timestamp.

### Reactome — pathway membership
- **Source**: Reactome ContentService
- **Citation**: Milacic M. *et al.* (2024) *Nucleic Acids Research*
  52(D1):D672-D678. doi:10.1093/nar/gkad1025

### HGNC — gene-symbol authority
- **Source**: https://rest.genenames.org/
- **Citation**: Seal R.L. *et al.* (2023) *Nucleic Acids Research*
  51(D1):D1003-D1009. doi:10.1093/nar/gkac888

## Cache layout

All downloads land under `$XDG_CACHE_HOME/pd_target_credentialing/`
(or `~/.cache/pd_target_credentialing/` if `XDG_CACHE_HOME` is unset).

## Provenance

Every result PNG in `results/` ships with a `.caption.md` sidecar recording
(a) the dataset version used, (b) the notebook cell that produced it,
(c) the commit SHA at which it was produced (§2.6).