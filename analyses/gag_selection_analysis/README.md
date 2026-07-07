# GAG biosynthesis selection analysis (lineage-specific selection scan)

This folder implements the **lineage-specific selection test** described in the
paper's Methods (subsection *"GAG biosynthesis selection analysis"*). The test
asks whether a group of regulatory elements (cCREs) associated with a given unit
of analysis — a single gene, or a functional gene set (e.g. a KEGG/GO/HPO term) —
shows a **coordinated directional regulatory shift** between the human-derived and
great-ape ancestral alleles that is more extreme than expected under a random model.

The method treats the per-cCRE MPRA log fold-changes (human-derived vs ancestral)
of a unit as a **cumulative trajectory** (a "random walk"). The endpoint of the
observed trajectory is compared against an empirical null distribution built by
repeatedly sampling random cCRE sets of matched size from the genome-wide pool.
An empirical *P*-value is the fraction of random trajectories whose endpoint is at
least as extreme as the observed one, and *P*-values are FDR-corrected
(Benjamini–Hochberg) across all tested units.

> Note: the engine is general-purpose. The pool-building notebook provided here
> builds pools at the **gene level**; the same `02`/`03` scripts run unchanged on
> pools built for functional gene sets (KEGG/GO/HPO) by supplying the appropriate
> `--id_col`.

---

## Pipeline & files

Run the three steps in order:

| Step | File | Role |
| :--- | :--- | :--- |
| 1 | `01_create_pools_gene_level.ipynb` | Build the experimental (`pool_genes_*`) and background (`pool_random_*`) pools from the master MPRA annotation table. |
| 2 | `02_random_trajectories_analysis.py` | Main engine: generate the random-trajectory null matrix, compute empirical *P*-values, FDR-correct, and produce diagnostic plots. |
| 3 | `03_replot_viz.py` | High-memory utility that re-plots publication-quality trajectory figures from the cached null matrix (no re-simulation). |
| 4 | `04_gag_pathway_selection_fdr.ipynb` | FDR-corrects the selection-test *P*-values across the four KEGG GAG pathways (added in revision; Reviewer #4 comment 6). |
| — | `plot_utils.py` | Shared plotting helpers used by the notebook. |

### Step 1 — Build the pools (`01_create_pools_gene_level.ipynb`)
Set the configuration at the top of the notebook and run all cells:
* `ANNOTATIONS_CSV` — path to the master per-oligo MPRA annotation table
  (`humanMPRA_annotations_v3.csv`).
* `min_n` — retain only genes linked to at least this many differentially active cCREs (integer ≥ 0).
* `condrocytes_filter` — restrict to cCREs accessible in chondrocyte ATAC-seq (`True`/`False`).
* `top_linked_filter` — restrict to high-confidence cCRE–gene links, i.e. ≥2 linking methods or within a promoter (`True`/`False`).

The notebook filters oligos (≥51 DNA counts in both alleles), explodes cCREs to a
single-gene-per-row table, **sets logFC to 0 for non-differentially-active cCREs**,
and writes two CSVs into a run folder named by the configuration
(`Run_minN_{X}_CondroFilt_{True/False}_TopLinkFilt_{True/False}/`):

* **`pool_genes_*.csv`** — one row per gene. Key columns:
  * `NCBI_Gene_ID`
  * `logFC_active` — list of per-cCRE logFC for the gene's active cCREs (non-diff-active set to 0)
  * `logFC_active_scuffled` — the same list, randomly permuted
  * `logFC_active_scuffled_cum` — cumulative sum of the permuted list (the observed trajectory)
* **`pool_random_*.csv`** — the genome-wide background pool. Key columns:
  * `oligo_ID`
  * `logFC_derived_vs_ancestral`

### Steps 2 & 3 — Run the scan and visualize
See the execution templates below. Step 2 caches the simulated null as a NumPy
matrix (`matrix_traj_*.npy`); step 3 reuses that cache. These `.npy` files are
large (often >10 GB) and are **git-ignored** (see `.gitignore`) — they are a local
cache, not a version-controlled artifact.

---

## Execution guide

The commands below are **generic templates** (LSF `bsub`); adjust queue, memory,
and statistical flags to your dataset. Memory scales with the number of
permutations and the maximum trajectory length: 32–64 GB is usually enough for
100K iterations; 128 GB+ is recommended for 500K iterations on large gene sets.

### Step 2 — Main analysis
```bash
bsub -q <QUEUE_NAME> -R "rusage[mem=<MEM_MB>]" -M <MEM_MB> \
  "python 02_random_trajectories_analysis.py \
  --genes_file <PATH_TO_POOL_GENES_CSV> \
  --pool_file  <PATH_TO_POOL_RANDOM_CSV> \
  --output_dir <OUTPUT_FOLDER_NAME> \
  --num_random <NUMBER_OF_PERMUTATIONS> \
  --pool_col logFC_derived_vs_ancestral \
  --id_col   <ID_COLUMN_NAME> \
  [--with_replacement] \
  [--min_nonzero <MIN_N>]"
```
Use a large `--num_random` for real runs (e.g. 100,000–500,000; the script default
of 1,000 is only for quick tests). For gene-set analyses (KEGG/GO/HPO), set
`--id_col` to the term-ID column and use `--with_replacement` (one cCRE may target
multiple genes within a term).

### Step 3 — Enhanced visualization
```bash
bsub -q short -R "rusage[mem=<MEM_MB>]" -M <MEM_MB> \
  "python 03_replot_viz.py \
  --genes_file   <PATH_TO_POOL_GENES_CSV> \
  --results_file <PATH_TO_analysis_results.csv> \
  --matrix_file  <PATH_TO_CACHED_NPY_MATRIX> \
  --output_plot  <PATH_TO_SAVE_IMAGE.png> \
  --id_col       <ID_COLUMN_NAME> \
  --sig_col      'reject_fdr_0.05'"
```

### Key parameters
| Parameter | Description |
| :--- | :--- |
| `--with_replacement` | Sample the background with replacement. Recommended for gene-set (KEGG/GO/HPO) analyses, where one cCRE can target several genes in a term. |
| `--min_nonzero` | Require each random trajectory to contain at least this many differentially active (non-zero) cCREs, matching the experimental signal density. |
| `--num_random` | Number of permutations building the null. Larger = more precise small *P*-values (and more memory/time). |
| `--id_col` | Column identifying the unit of analysis (`NCBI_Gene_ID` for gene level; a term ID for gene sets). |

---

## Outputs (per run directory)
* **`analysis_results.csv`** — per-unit table: empirical `p_value`, FDR-adjusted `p_adj`, `reject_fdr_0.05`, trajectory length `L`, and endpoint.
* **`trajectories_plot.png`** — observed trajectories overlaid on the null "cloud".
* **`qq_plot_log.png`** — −log10 QQ plot with a Kolmogorov–Smirnov test of the *P*-value distribution.
* **`pvalue_histogram.png`** — raw *P*-value histogram vs the uniform expectation.
* **`matrix_traj_*.npy`** — cached null matrix (git-ignored; reused by step 3).

---

## Requirements
Python 3.x with: NumPy, pandas, Matplotlib, statsmodels, SciPy (see
`requirements.txt`). The pool-building notebook additionally imports `plot_utils.py`
from this folder.

---

