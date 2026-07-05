# Figure generation

Each notebook produces the raw panels for one figure of the paper (assembled/labelled
afterwards in a vector editor). The notebooks are intentionally **standalone** — each
loads and preprocesses the MPRA data it needs, so it can be run on its own. This means
some loading/normalization code is repeated across the QC notebooks by design.

| Notebook | Produces |
| :--- | :--- |
| `produce_hMPRA_fig1.ipynb` | Figure 1 (library design & MPRA QC), incl. the activity-vs-open-chromatin ECDF |
| `produce_hMPRA_fig3.ipynb` | Figure 3 (skeletal ECM / GAG): functional-enrichment dot plots (MPRA + hybrid ASE) and the GAG pathway panel |
| `produce_hMPRA_fig_DMMB.ipynb` | Figure 4 (DMMB GAG assay): per-joint box plots, human-vs-ape fold changes, age/sex normalization and validations |
| `produce_hMPRA_fig5.ipynb` | Figure 5 (*EVC* / Hedgehog): *EVC* barcode violin/swarm/box plots and the genome-wide 5-Mb window enrichment scan |
| `produce_hMPRA_supp_fig_MPRA_qc.ipynb` | Extended Data QC — **differential**-activity: volcano, differential-activity density, control scatters, GC-content bias, cross-library positive-control correlation |
| `produce_hMPRA_supp_fig_MPRA_qc_2.ipynb` | Extended Data QC — **activity**: replicate correlation, RNA-vs-DNA, activity distribution, ancestral-vs-derived allele activity |
| `produce_hMPRA_stats_for_paper.ipynb` | Numeric values quoted in the text (counts of active / differentially active cCREs, controls, library totals, top-candidate details); produces no figures |

## Dependencies & configuration

* **Shared helpers:** all notebooks import `const` (figure styling, `save_fig`,
  colors) from `../common/const.py`; each adds `analyses/common` to `sys.path`.
* **Paths (edit at the top of each notebook):**
  * `output_path` — directory where the panels are written (PNG + SVG).
  * `WORK_DIR` — working directory used for the few relative-path reads/writes
    (`fig5` reads `active.bed` / `diff_active.bed`; `fig3` writes `hybrid_data_*.csv`).
* **Inputs:** the notebooks read the MPRA activity / comparative tables, the master
  annotation table `humanMPRA_annotations_v3.csv`, the hybrid `ExpLBM_polarization_results.tsv`,
  enrichment outputs, SCREEN cCRE beds, and the DMMB assay table — all from the lab
  data tree (absolute paths).
* Notebook outputs were cleared before committing; re-run to regenerate the panels.
