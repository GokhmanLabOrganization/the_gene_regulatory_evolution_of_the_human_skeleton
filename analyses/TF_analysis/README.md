# Transcription factor binding analysis

This folder analyses how transcription factor (TF) binding differs between the
human-derived and great-ape ancestral alleles of the MPRA cCREs, and links that
differential binding to differential MPRA activity (Methods, *"Transcription
factor binding sites"*). Two orthogonal approaches are used: **FIMO** (predicted
motif matches, JASPAR) and **PBM** (protein binding microarray affinities). The
analysis identifies TFs that are both enriched in differentially active cCREs and
whose differential binding correlates with differential activity.

## Scope

The **differential-binding scoring itself** (the FIMO/PBM computational engine) is
a separate tool and is **not** part of this repository — see
[GokhmanLabOrganization/differential-TF-binding](https://github.com/GokhmanLabOrganization/differential-TF-binding).
This folder contains the project-specific **input preparation** for that tool and
the **downstream analysis** of its output.

## Pipeline

Input preparation (feeds the external tool):
| File | Produces |
| :--- | :--- |
| `create_FIMO_input.ipynb` | variant BED file for the tool |
| `filter_jaspar_tfs.ipynb` | JASPAR MEME motif file filtered to TFs expressed in chondrocytes (for FIMO) |
| `filter_PBM_input.ipynb` | PBM 8-mer data filtered to chondrocyte-expressed TFs |
| `config_file_chonds.tsv` | configuration for the differential-TF-binding tool |

↓ **run the external `differential-TF-binding` tool** → `{FIMO,PBM}_diff_binding_scores.tsv`

Analysis of the tool output:
| File | Role |
| :--- | :--- |
| `process_fimo_tool_output.ipynb` | merge the tool's differential-binding scores with the MPRA oligo annotations |
| `analyze_processed_fimo_tool_output.ipynb` | main analysis: motif enrichment in differentially active cCREs, correlation of differential binding with differential activity (FIMO + PBM), and the set of TFs significant across tests |
| `analyze_modern_vs_archaic.ipynb` | compare TF binding signals between this (great-ape) MPRA and the modern-human MPRA (osteoblasts / NPCs) |
| `human_chimp_trans_env_analysis.ipynb` | *trans*-environment validation (e.g. concordance of shared positive controls across MPRAs) |

## Dependencies & configuration

* **Shared helpers:** the analysis notebooks import `const` (figure styling,
  `save_fig`) from `../common/const.py`; each adds `analyses/common` to `sys.path`.
* **Working directory:** each notebook sets `BASE_DIR` at the top (the lab
  `Collaboration` data root) and all relative paths resolve from there.
* **Inputs:** the master MPRA table `humanMPRA_annotations_v3.csv`, JASPAR 2024
  motifs, PBM 8-mer data, and cell-type expression tables (chondrocyte / osteoblast
  / NPC), all under the lab data tree.
* Notebook outputs were cleared before committing; re-run to regenerate figures.
