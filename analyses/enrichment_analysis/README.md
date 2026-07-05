# Functional enrichment analyses

This folder runs the paper's **functional enrichment analyses** (Methods,
*"Functional enrichment analyses"*): testing whether genes linked to regulatory
changes are enriched for particular biological functions, phenotypes, and diseases,
using GO, KEGG, the Human Phenotype Ontology (HPO), the Genetic Association Database
(GAD), and Gene ORGANizer.

Two entry-point R Markdown wrappers, one per gene set of interest:

| File | Gene list tested | Background |
| :--- | :--- | :--- |
| `01_enrich_mpra_diff_active.Rmd` | Genes linked (via elite associations) to **MPRA differentially active** cCREs | Genes linked to active cCREs |
| `02_enrich_hybrid_ase.Rmd` | Genes with **allele-specific expression (ASE)** in the human-ape hybrids | Genes tested for ASE (ASE + nonASE) |

Each wrapper runs an enrichment test (gene list vs background) across all databases,
and `01` additionally runs a "compare" test (up- vs down-regulated genes).

## The enrichment toolkit

The wrappers call the lab's enrichment toolkit, a **vendored copy** of
[GokhmanLabOrganization/Enrichment-analysis](https://github.com/GokhmanLabOrganization/Enrichment-analysis),
bundled here under `enrichment_toolkit/`. Each wrapper sources it via
`grandEnrichDir <- "enrichment_toolkit"`. The toolkit provides
`TheGrandEnrich_*` (GO/KEGG/GAD/…), `HPO_*`, and `ORGANizer_*` functions; see
`enrichment_toolkit/README.md` for its own documentation, arguments, and
`environment.yml`.

> Note: the toolkit reads the underlying ontology/database files from the lab
> environment (`locORser = "ser"` for server, `"loc"` for local). Those database
> files are not part of this repository.

## Inputs (set at the top of each wrapper)

`01_enrich_mpra_diff_active.Rmd`:
* `data_dir` — master MPRA per-oligo annotation table (`humanMPRA_annotations_v3.csv`).
* `base_dir` — where run outputs (`base_dir/runs/<run_id>/`) and the master log are written.

`02_enrich_hybrid_ase.Rmd`:
* `data_dir` — per-gene ASE summary from the hybrids (`ASE_info.tsv`).
* `ASE_CELLTYPE` — which cell type's ASE genes to test (the table has one
  `*_gene_ase_type` column per cell type). **Confirm this matches the paper's
  osteochondral progenitor / chondrocyte column.**
* `base_dir` — output location, as above.

Both wrappers restrict to elite cCRE-gene associations and to genes expressed in
the relevant cell type, matching the Methods.

## Outputs (per run directory)
* `<run_id>_final_results.csv` — combined enrichment table (Term, FoldChange, P-value, FDR) across databases, filtered to FDR < 0.05.
* `<run_id>_gene_list.txt/.csv`, `<run_id>_background.txt/.csv` — the input gene sets.
* `<run_id>_parameters.txt` — the parameters used for the run.
* `*_enrichment_results*.png` — bar-plot visualizations.
* A master log (`enrich_master_table.csv`) accumulating one row per run.

## Requirements
R with `dplyr`, `tidyr`, `ggplot2`, `stringr` (wrappers) plus the toolkit's own
dependencies (see `enrichment_toolkit/environment.yml`; `readxl`).
