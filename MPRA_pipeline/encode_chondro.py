# todo: think about how to treat spike ins - remove! but when?
# remove number after dot in ensembl id

import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import venn
import numpy as np
from itertools import chain
import seaborn as sns

# load datasets
primary_1 = pd.read_csv("./region_gene_link/chondrocytes/expression_chondrocytes/encode/ENCSR000CUE/ENCFF023BMB.tsv", sep='\t', header=0)
primary_2 = pd.read_csv("./region_gene_link/chondrocytes/expression_chondrocytes/encode/ENCSR000CUE/ENCFF633DRY.tsv", sep='\t', header=0)
invitro_1 = pd.read_csv("./region_gene_link/chondrocytes/expression_chondrocytes/encode/ENCSR774MGO/ENCFF054ADN.tsv", sep='\t', header=0)
invitro_2 = pd.read_csv("./region_gene_link/chondrocytes/expression_chondrocytes/encode/ENCSR774MGO/ENCFF217MKR.tsv", sep='\t', header=0)
invitro_3 = pd.read_csv("./region_gene_link/chondrocytes/expression_chondrocytes/encode/ENCSR774MGO/ENCFF496AMJ.tsv", sep='\t', header=0)

expressed_genes = []
plt.clf()
fig, axes = plt.subplots(5, figsize=(4, 20), sharex=True)
for n, dataset in enumerate([primary_1, primary_2, invitro_1, invitro_2, invitro_3]):
    #filter our spikeins
    dataset = dataset[~dataset['gene_id'].str.contains('Spike')]
    
    # draw hist
    sns.histplot(ax=axes[n], data=dataset, x = "TPM", bins=30)
    
    # filter for tmp >1 and count
    expressed = dataset[(dataset['TPM'] > 1)] 
    genes = set(expressed["gene_id"].tolist()) #turn into list using only unique oligos
    print(len(genes))
    expressed_genes.append(genes)

plt.suptitle(f'Histogramm of TMP for different datasets')
plt.savefig(f'./region_gene_link/chondrocytes/expression_chondrocytes/encode/hist_tpm.pdf')
plt.savefig(f'./region_gene_link/chondrocytes/expression_chondrocytes/encode/hist_tpm.png', dpi=330)

# draw venn diagramm 
labels = venn.get_labels(expressed_genes, fill=['number', 'logic'])
plt.clf()
fig, ax = venn.venn5(labels, names=['primary_1', 'primary_2', 'invitro_1', 'invitro_2', 'invitro_3'])
plt.title(f'Expressed genes - overlap')
plt.savefig(f'./region_gene_link/chondrocytes/expression_chondrocytes/encode/venn_expressed.pdf')
plt.savefig(f'./region_gene_link/chondrocytes/expression_chondrocytes/encode/venn_expressed.png', dpi=330)

# combine all expressed genes, record number
expressed_genes_combined = list(set(chain(*expressed_genes)))
print(f"expressed genes in combined datasets: {len(expressed_genes_combined)}")

# remove the point from ensemble ids
processed_ids = []
for id_unprocessed in expressed_genes_combined:
    id_processed = id_unprocessed.split('.')[0]
    processed_ids.append(id_processed)
    
# save expressed genes in format suitable for genealacart
with open(r'./region_gene_link/chondrocytes/expression_chondrocytes/encode/expressed_genes_mixed_ids.txt', 'w') as file:
    for item in processed_ids:
        file.write("%s\n" % item)


# out side of script - use genealacart to convert to ncbi ids

# other scirpt - do set on ncbi ids, save ncbi ids in format suitable for yardens script