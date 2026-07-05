# Functional enrichment analyses

This folder runs the paper's **functional enrichment analyses** (Methods,
*"Functional enrichment analyses"*): testing whether genes linked to regulatory
changes are enriched for particular biological functions, phenotypes, and diseases,
using GO (Biological Process), KEGG, the Genetic Association Database (GAD), the
Human Phenotype Ontology (HPO), and Gene ORGANizer.

Two entry-point R Markdown wrappers, one per data source:

| File | Gene list tested | Background |
| :--- | :--- | :--- |
| `01_enrich_mpra_diff_active.Rmd` | Genes linked to **MPRA differentially active** cCREs (`\|FC\| > 1.333`, elite associations) | Genes linked to active cCREs (≥50 DNA counts in both alleles) |
| `02_enrich_hybrid_ase.Rmd` | **Human-derived** ASE genes in the hybrids (ASE ∩ `derived == 'human-derived'`) | Genes tested for ASE (ASE + nonASE) |

Each wrapper runs two tests:
* **diff_vs_active** — the gene list vs the background, across all databases.
* **compare** — up- vs down-regulated genes against each other (relative enrichment).

## The enrichment toolkit

The wrappers call the lab's enrichment toolkit, a **vendored copy** of
[GokhmanLabOrganization/Enrichment-analysis](https://github.com/GokhmanLabOrganization/Enrichment-analysis),
bundled here under `enrichment_toolkit/`. Each wrapper sources it via
`grandEnrichDir <- "enrichment_toolkit"`. It provides the `TheGrandEnrich_*`
(GO/KEGG/GAD), `HPO_*`, and `ORGANizer_*` enrich/compare functions; see
`enrichment_toolkit/README.md` for its own docs and `environment.yml`.

> Note: the toolkit reads the underlying ontology/database files from the lab
> environment (`locORser = "ser"` for server, `"loc"` for local). Those database
> files are not part of this repository.

## Inputs (set at the top of each wrapper)

`01_enrich_mpra_diff_active.Rmd`:
* `data_dir` — master MPRA per-oligo annotation table (`humanMPRA_annotations_v3.csv`).
* `base_dir` — output root; results go under `base_dir/<parameters$name>/`.
* `parameters$diff_threshold` (1.333) and `parameters$compare_threshold` (1.333) — fold-change cutoffs (applied as `log2()`); `min_dna_counts_*` = 50.

`02_enrich_hybrid_ase.Rmd`:
* `data_dir` — hybrid ASE + lineage-polarization table (`ExpLBM_polarization_results.tsv`; `ExpLBM` = limb-bud mesenchyme / osteochondral progenitors). Uses the `ExpLBM_gene_ase_type` (ASE/nonASE) and `derived` (human-/chimp-derived) columns.
* `base_dir` — output root, as above.
* `parameters$diff_threshold` (0.5), `compare_threshold` (0.2), `min_TPM` (0).

Both wrappers create their output subdirectories automatically (`dir.create`).

## Outputs (under `base_dir/<parameters$name>/`)
* `diff_vs_active/` — per-database raw enrichment tables (`kegg_raw_output.csv`,
  `GAD_disease_raw_output.csv`, `GO_raw_output.csv`, `HPO_raw_output.csv`,
  `ORGANizer_raw_output.csv`), a combined `diff_vs_active_processed_results.csv`
  (FDR < 0.05), and the gene/background lists.
* `compare/` — the equivalent raw + `compare_processed_results.csv` for the
  up-vs-down comparison, plus the up/down gene lists.
* `parameters.csv` — the parameters used for the run.

## Requirements
R with `dplyr`, `tidyr`, `ggplot2`, `stringr` (wrappers) plus the toolkit's own
dependencies (see `enrichment_toolkit/environment.yml`; `readxl`).
