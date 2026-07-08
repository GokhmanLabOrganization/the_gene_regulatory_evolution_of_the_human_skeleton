# The Gene Regulatory Evolution of the Human Skeleton

Code accompanying the manuscript *"The Gene Regulatory Evolution of the Human
Skeleton"* (Gokhman Lab, Weizmann Institute of Science; Inoue Lab, Kyoto
University). The paper has been accepted but is **not yet published**; this
repository will be finalized upon publication.

## Overview

The study maps the gene-regulatory changes that emerged during human skeletal
evolution. Massively parallel reporter assays (MPRAs) in chondrocytes assayed
561,410 fixed or nearly-fixed human-derived substitutions in promoters and
enhancers, and human–ape hybrid cells differentiated into osteochondral
progenitors were used to separate *cis*- from *trans*-regulatory effects. Together
they yield a genome-wide atlas of human-specific *cis*-regulatory changes,
highlighting divergence in Hedgehog signaling and extracellular-matrix /
glycosaminoglycan (GAG) pathways.

## Repository structure

```
.
├── MPRA_pipeline/              # Processing of the MPRA data (DNA/RNA counts -> regulatory activity)
├── human-ape_hybrid_pipeline/  # Allele-specific expression in human–ape hybrid cells (cis/trans)
└── analyses/                   # Downstream analyses and paper figures
    ├── common/                 # Shared plotting / config helpers (const.py)
    ├── enrichment_analysis/    # Pathway / gene-set enrichment of MPRA-active loci and hybrid ASE
    ├── gag_selection_analysis/ # Lineage-specific selection scan on GAG / ECM pathways
    ├── TF_analysis/            # Differential transcription-factor binding (FIMO / PBM)
    ├── MPRA_LegNet/            # Sequence-to-activity model (LegNet) predictions
    ├── figure_generation/      # Notebooks producing the main and supplementary figures
    └── additional_analyses/    # Supporting analyses added during peer review
```

Each folder under `analyses/` contains its own README describing its inputs,
outputs, and scripts.

## Reproducing the analyses

Input-data paths are surfaced at the top of each script/notebook (a
`# --- Paths (EDIT) ---` config block, or a "set paths here" section) and must be
edited to point at local copies of the data. Raw and processed data will be made
available upon publication; see the paper's Data Availability statement.
