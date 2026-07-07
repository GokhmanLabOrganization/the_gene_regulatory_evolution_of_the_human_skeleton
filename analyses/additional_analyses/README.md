# Additional analyses (revision)

Analyses added during peer review, each answering a specific reviewer comment.
Every notebook begins with a markdown cell stating the comment it addresses and a
short summary of our response. The notebooks are standalone and share the plotting
helpers in `../common/const.py`; two of them also use an external liftOver tool
(surfaced as `LIFTOVER_DIR` at the top of the notebook).

| Notebook | Reviewer comment | What it does |
| :--- | :--- | :--- |
| `variant_fixation_1KG_HGDP.ipynb` | R3 minor + R4 #1 | Re-check allele frequencies of the assayed variants in gnomAD v4 and 1KG+HGDP (v3.1) |
| `barcode_outlier_filtering_stats.ipynb` | R3 minor | Quantify how many barcodes the extreme-RNA-outlier filter removes |
| `seq325030_subsampling.ipynb` | R3 minor | Document the barcode subsampling for seq_325030 |
| `active_cCRE_open_chromatin_fisher.ipynb` | R4 #3 | Exact active/non-active × open/closed-chromatin counts + Fisher's exact test |
| `cCRE_atacseq_cell_type_specificity.ipynb` | R3 minor | Correlate cCRE activity with chromatin accessibility across cell types (Ext. Data Fig. 3) |
| `hybrid_expression_checks.ipynb` | R3 | Chondrocyte-vs-osteochondral expression correlation (r=0.745) and gorilla-specific gene count |
| `trans_environment.ipynb` | R3 | Compare TF expression between human and chimpanzee (*trans*-environment similarity) |
| `osteoarthritis_kegg_enrichment.ipynb` | R4 (within-human variation & disease) | KEGG-pathway enrichment of osteoarthritis-associated genes |
| `eqtl_gwas_overlap.ipynb` | R4 (within-human variation & disease) | eQTL and GWAS-catalog overlap of cCREs, stratified by MPRA activity |

A related revision analysis — the multiple-testing correction of the lineage-specific
selection test on the four KEGG GAG pathways (Reviewer #4 comment 6) — lives with the
selection scan it belongs to: `../gag_selection_analysis/04_gag_pathway_selection_fdr.ipynb`.

## Provenance

These notebooks were split out of four combined "miscellaneous" notebooks so that each
analysis stands on its own and maps to a single reviewer comment. Notebook outputs were
cleared before committing; re-run to regenerate figures/tables. Inputs are read from the
lab data tree (absolute paths); set `output_path` (and `LIFTOVER_DIR` where present) at the
top of each notebook.
