# Evolutionary Regulatory Trajectory Analysis (MPRA)

This project is designed to identify coordinated directional selection in human regulatory elements using Massively Parallel Reporter Assay (MPRA) data. The core methodology employs a **Stochastic Trajectory Model** to determine if a group of enhancers associated with a specific functional category (Genes,KEGG, HPO, or GO) shows a significant cumulative regulatory shift compared to an empirical null distribution.

---
## 📂 Project Directory Structure

The repository is organized into hierarchical levels based on the unit of analysis (Genes vs. Pathways) and the specific functional database used.

### 1. Analysis Hierarchy
* **`gene_level/`**: Simulations where the stochastic trajectory is calculated for individual genes based on their linked enhancers.
* **`pathways_level/`**: Aggregated analyses where trajectories represent functional groups. This is subdivided into:
    * `kegg/`: Kyoto Encyclopedia of Genes and Genomes.
    * `GO-BP` / `GO-CC`: Gene Ontology (Biological Process / Cellular Component).
    * `hpo_genes_to_phenotype/`: Human Phenotype Ontology (Skeletal-focused curation).

### 2. Standardized Execution Folders
Each simulation run is stored in a dedicated directory with a naming convention that reflects its hyper-parameters:
`Run_minN_{X}_CondroFilt_{True/False}_TopLinkFilt_{True/False}_WithReplace_{True/False}`

* **`minN`**: Minimum number of differentially active oligos required per pathway.
* **`CondroFilt`**: Whether tissue-specific (Chondrocyte) filtering was applied.
* **`TopLinkFilt`**: Whether high-confidence genomic linkage criteria were enforced.
* **`WithReplace`**: The sampling strategy used for the background null model.

### 3. Standardized Output Files
Every results directory contains the following core components:
* **`analysis_results.csv`**: The primary data table containing Pathway IDs, empirical P-values, and Benjamini-Hochberg adjusted P-values (FDR).
* **`trajectories_plot.png`**: A visualization of the experimental pathways (purple lines) compared against a subset of the null distribution (light blue cloud).
* **`qq_plot_log.png`**: A -log10 QQ-plot with Kolmogorov-Smirnov (KS) test results to validate P-value distribution.
* **`pvalue_histogram.png`**: A frequency distribution of raw P-values against the expected uniform null.
* **`pool_genes_*.csv` / `pool_random_*.csv`**: The input data pools generated for that specific configuration.

### ⚡ Optimization & Heavy Files
* **`.npy` Matrices**: The main simulation script (`random_trajectories_analysis.py`) generates a serialized NumPy matrix of 500,000 random trajectories (`matrix_traj_*.npy`). 
* **Note**: These files are extremely heavy (>10GB) and are typically excluded from version control. They serve as a local cache to allow `replot_viz.py` to generate high-resolution plots instantaneously without re-running the full simulation.

### 📜 Core Scripts
* `random_trajectories_analysis.py`: The main engine for stochastic simulations and vectorized P-value computation.
* `replot_viz.py`: A high-memory utility for generating publication-quality figures from cached matrices.
* `plot_utils.py`: Shared plotting library for consistent visual styles across the project.
* `create_pools_*.ipynb`: Data preparation notebooks for database merging and link expansion ("exploding").

## 🛠 Project Workflow

The analysis is divided into two main phases:

### 1. Data Preparation (`create_pools_*.py`)
Dedicated scripts for each functional database (GENES, KEGG, HPO, GO CC, GO BP) process the raw MPRA data. 
* **Biological Filtering:** Restricts enhancers to those show availability in ATAC seq in chondrocytes and applies high-confidence gene-linkage filters.
* **LogFC Neutralization:** Sets $logFC = 0$ for non-differentially active oligos to focus the signal on high-confidence evolutionary divergence.
* **Functional Mapping:** Performs an inner join between the filtered MPRA library and the chosen functional database.
* **Output:** Generates `pool_genes_*.csv` (experimental trajectories) and `pool_random_*.csv` (the background sampling universe).



### 2. Statistical Simulation (`random_trajectories_analysis.py`)
This is the main analysis engine. It utilizes a high-performance vectorized approach to run large-scale permutations.
* **Null Model Generation:** Constructs a matrix of hundreds of thousands of random trajectories sampled from the background pool.
* **Empirical P-value Calculation:** Compares the final cumulative sum of the experimental trajectory against the random distribution at the same step length ($L$).
* **FDR Correction:** Applies the Benjamini-Hochberg procedure to control for false discoveries across all tested categories.

### 3. Visualization Utility: `replot_viz.py`
This script is used for high-quality, publication-ready visualization of the results. Instead of recalculating the null model, it leverages the cached `.npy` matrix to plot statistical thresholds (P=0.05) and differentiate significant vs. non-significant trajectories with high visual contrast.

#### **Key Features**
* **Statistical Accuracy:** Calculates 95th percentile thresholds using the *full* 500k matrix.
* **Visual Clarity:** Plots only a representative subset of the background (10,0000 trajectories) to prevent "color bloating" while maintaining context.
* **Significance Highlighting:** Automatically colors trajectories based on FDR results (Navy for significant, DodgerBlue for non-significant).

---

## 🔬 The Stochastic Trajectory Model

For each functional category, we calculate the cumulative sum of the human allelic effect ($logFC$) compared to the ancestral state:

$$S_n = \sum_{i=1}^{n} \log_2\left(\frac{\text{Derived}}{\text{Ancestral}}\right)_i$$

A significant deviation from the random "cloud" indicates a coordinated gain or loss of regulatory activity in the human lineage for that specific biological program.


---

## 💻 Execution Guide 

### Phase 1: Main Analysis
The following commands are **generic templates**. You should adjust the memory allocation and statistical flags based on your specific research question and dataset size.

### Generic BSUB Template
```bash
bsub -q <QUEUE_NAME> -R "rusage[mem=<MEM_MB>]" -M <MEM_MB> \
  "python random_trajectories_analysis.py \
  --genes_file <PATH_TO_GENES_POOL_CSV> \
  --pool_file <PATH_TO_RANDOM_POOL_CSV> \
  --output_dir <OUTPUT_FOLDER_NAME> \
  --num_random <NUMBER_OF_PERMUTATIONS> \
  --id_col <ID_COLUMN_NAME> \
  [--with_replacement] \
  [--min_nonzero <MIN_N>]"

### Phase 2: Enhanced Visualization (`replot_viz.py`)
Use this generic template to generate high-resolution plots after the main analysis is complete. Note that this script requires high memory (32GB-150GB) depending on the size of the `.npy` matrix.

```
bsub -q short -R "rusage[mem=<MEM_MB>]" -M <MEM_MB> \
  "python replot_viz.py \
  --genes_file <PATH_TO_POOLS_GENES_CSV> \
  --results_file <PATH_TO_ANALYSIS_RESULTS_CSV> \
  --matrix_file <PATH_TO_CACHED_NPY_MATRIX> \
  --output_plot <PATH_TO_SAVE_IMAGE.png> \
  --id_col <ID_COLUMN_NAME> \
  --sig_col 'reject_fdr_0.05'"

  for example: 
  bsub -q short -R "rusage[mem=64000]" -M 64000 \
  "python /home/labs/davidgo/idanko/backup/VScode/script/random_trajectories_analysis.py \
  --genes_file ./pool_genes_HPO_minN_0_CondroFilt_False_TopLinkFilt_False_WithReplace_False.csv \
  --pool_file ./pool_random_HPO_minN_0_CondroFilt_False_TopLinkFilt_False_WithReplace_False.csv \
  --output_dir results_hpo \
  --num_random 100000 \
  --pool_col logFC_derived_vs_ancestral \
  --id_col hpo_id \
  --with_replacement \
  --min_nonzero 1" 

### Optional Parameters & Logic

The analysis script provides several flags to adjust the statistical model according to the biological context:

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `--with_replacement` | **Optional** | **Sampling with replacement:** If enabled, a single enhancer can be sampled multiple times for a random trajectory. This is highly recommended for GO/KEGG/HPO analyses where one enhancer might target multiple genes within the same functional category (capturing "Regulatory Dosage"). |
| `--min_nonzero` | **Optional** | **Minimum Diff-Active Threshold:** Ensures every random trajectory contains at least $X$ differentially active enhancers. This ensures the Null Model matches the experimental pathways' complexity and signal density. |
| **Memory Allocation** | **Generic** | Memory usage scales with the number of iterations and max trajectory length. For 100K iterations, 32GB-64GB is usually sufficient. For 500K iterations on large Gene Ontology terms, 128GB+ is recommended. |

---

## 📊 Outputs & Diagnostics

The script generates a results directory for each run containing the following files:

### 1. Statistical Results
* **`analysis_results.csv`**: The primary output table. It includes raw P-values, FDR-corrected P-values (`p_adj`), the number of enhancers ($L$), and the final cumulative endpoint for each category.

### 2. Visualization
* **`trajectories_plot.png`**: A visualization of the "Random Walks." The experimental trajectories (purple) are overlaid on the background null distribution "cloud" (blue).


* **`qq_plot_log.png`**: A -log10 QQ-plot used to diagnose P-value distribution. Significant deviations from the diagonal line indicate biological signals beyond random expectation. It also includes the **Kolmogorov-Smirnov (KS) test** results.


* **`pvalue_histogram.png`**: A distribution plot of raw P-values to check for enrichment at the lower end ($P < 0.05$).

### 3. Optimization Cache
* **`matrix_traj_*.npy`**: A serialized NumPy matrix of the random permutations. If you re-run the analysis with the same parameters, the script will load this cache, making the re-analysis near-instantaneous.

---

## ⚙️ Technical Requirements

The pipeline requires a Python 3.x environment with the following libraries:

* **NumPy & Pandas**: For high-performance data manipulation and matrix operations.
* **Matplotlib & Seaborn**: For generating diagnostic plots.
* **Statsmodels**: Specifically for `multipletests` (FDR correction) and `qqplot`.
* **SciPy**: For performing the Kolmogorov-Smirnov (KS) test.
* **ast (literal_eval)**: For parsing list structures stored in CSV files.

---

## 🔬 Scientific Context
This pipeline was developed to analyze the evolutionary history of human-specific regulatory changes. By treating regulatory input as a cumulative trajectory, we can detect selection on biological *processes* even when individual enhancers show modest effects.

**Contact:** Idan Korenfeld | Gokhman Lab | Weizmann Institute of Science.