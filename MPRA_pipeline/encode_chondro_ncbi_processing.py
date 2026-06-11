import pandas as pd

ncbi_ids = pd.read_csv("./region_gene_link/chondrocytes/expression_chondrocytes/encode/genealacart_results_combined.csv", header=0)
print(ncbi_ids.shape[0])
ncbi_ids.dropna(inplace=True)
print(ncbi_ids.shape[0])
genes = set(ncbi_ids["NCBI Gene"].tolist())
print(len(genes))

# save expressed genes in format suitable for yardesn script
with open(r'./region_gene_link/chondrocytes/expression_chondrocytes/encode/expressed_genes_ncbi_ids.txt', 'w') as file:
    for item in genes:
        file.write("%s\n" % item)