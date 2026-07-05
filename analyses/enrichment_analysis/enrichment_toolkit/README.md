# 🧬 Enrichment Analysis

## 🚀 Usage

### 1. 📦 Clone the repository<b> 1. Clone the repository: </b>
   ```bash
   cd your_path
   git clone https://github.com/GokhmanLabOrganization/Enrichment-analysis
   cd Enrichment-analysis
   ```
<b> 2. 🔬 Enerichment scripts - Gene list vs. background: </b>
   ```These scripts test whether terms are enriched in a gene list of interest relative to a background gene set.
The gene list must be a subset of the background, or no background should be provided.

Available scripts:

TheGrandEnrich_enrich.R – GO, Reactome, GAD Disease, OMIM, and KEGG enrichment

HPO_enrich.R – Human Phenotype Ontology (HPO) enrichment

GeneORGANizer_enrich.R – Gene ORGANizer term enrichment
   ```

<b> 3. ⚖️ Compare nerichment scripts - Gene list 1 vs. Gene list 2: </b>
   ```
These scripts run enrichment analyses on two gene lists and then compare the enrichment values between them, calculating a relative enrichment score for each term.

This allows direct comparison of functional enrichment between two conditions or datasets.
   ```
<b>4. Output:</b>  
<pre><code>Each script produces a data frame with one row per term.
For each term, the output includes:

Number of genes associated with the term in each set

Expected number of genes

Fold change

P-value

Multiple-testing–corrected p-value (FDR), if the term passes the minimum gene threshold
</code></pre>

<b> 5. 🧩 Shared arguments: </b>
   ```
genelist - a gene set of interest 
background - a larger set of genes, must include all genes in genelist or empty
enrORdep - the direction of the test, can be enr for enrichment (higher probability to see term in genelist), dep (lower probability) or both (a non directional test) 
minGenes - the minimal number of gene associations to a term, if lower than the threshold, the script won't take the test into account (no FDR) 
Comp - a string that describes our test – can be “active vs. all” for example. In the compare scripts we have instead comp_main, comp1, and comp2
uniq - determines whether genelist will be tested with redundancy or not. If yes, a gene can appear multiple times and its associations will be counted separately
locORser - loc - running locally or ser - on the server
outdir/outpath - a path for saving the results file
   ```
<b> 6.📝  Other arguments: </b>
   ```
Arguments that are unique to specific scripts are documented within the corresponding R files.
   ```
<b> 7. ▶️ Example run: </b>
   ```
An example demonstrating how to run an enrichment script is provided in: run_example.R


   ```
<b> 8. 📚 Dependencies: </b>
```

- 🐍 **Conda**
  - Environment defined in `environment.yml`
  # Create the environment
  conda env create -f environment.yml

  # Activate the environment
  conda activate enrichment_env
  

- 📦 **R packages**
  - `readxl`

- 🧰 **Base R (no additional installation required)**

```

