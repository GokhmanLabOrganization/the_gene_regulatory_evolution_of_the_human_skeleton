# this script creates a summary table per gene to find top candaidates

import pandas as pd
import sys

cells = sys.argv[1]

mpra = pd.read_csv(f'./{cells}/comparative_analysis_combined/humanMPRA_with_seq_final2.csv', header=0, usecols=["oligo","RNA_DNA_ratio_log2_ancestral", "RNA_DNA_ratio_log2_derived", "DNA_counts_raw_ancestral","DNA_counts_raw_derived","barcode_count_ancestral", "barcode_count_derived", "normalized_activity_estimate_ancestral", "normalized_activity_estimate_derived", "activity_fdr_ancestral", "activity_fdr_derived", "activity_ancestral", "activity_derived", "logFC_derived_vs_ancestral", "differential_activity_fdr", "differential_activity"])
genes = pd.read_csv(f'./region_gene_link/{cells}/yardens_script/output_1/Summary_PerRegion2GeneAssociations.tsv', sep='\t', header=0)

# aggregate the mpra data per gene
mpra['id'] = mpra['oligo'].str.split('_').str[2]
genes['id'] = genes['chromosome'].astype(str) + ':' + (genes['start'] + 1).astype(str) + '-' + genes['end'].astype(str)
genes.dropna(subset=["Gene_symbol"], inplace=True)

merged_df = pd.merge(genes, mpra, on='id', how="inner")

def process_group(x):
    active_count_ancestral = (x['activity_ancestral'] == 'active').sum()
    active_count_derived =  (x['activity_derived'] == 'active').sum()
    diff_active_count = (x['differential_activity'] == True).sum()
    log_FC = x.loc[x['differential_activity'] == True, 'logFC_derived_vs_ancestral'].tolist()
    num_oligos = len(x)  
    return pd.Series({'num_oligos': num_oligos, 'active_count_ancestral': active_count_ancestral, 'active_count_derived': active_count_derived, 'diff_active_count': diff_active_count, 'log_FC': log_FC})
    
result = merged_df.groupby('Gene_symbol').apply(process_group).reset_index()
result_filtered = result[['Gene_symbol', 'num_oligos', 'active_count_ancestral', 'active_count_derived', 'diff_active_count', 'log_FC']]

print(result_filtered.head().to_string())

result_filtered['log_FC'] = result_filtered['log_FC'].apply(lambda x: sorted(x, reverse=True))

result_filtered['sum_of_logFC'] = result_filtered['log_FC'].apply(lambda x: sum(x))

max_length = max(len(lst) for lst in result_filtered['log_FC'])
for i in range(max_length):
    result_filtered[f'logFC_oligo{i+1}'] = result_filtered['log_FC'].apply(lambda x: x[i] if i < len(x) else None)
  
print(result_filtered.head().to_string())

# merge with ASE

ase = pd.read_csv(f'/home/labs/davidgo/Collaboration/data/Hybrids/ASE_info.tsv', header=0, sep='\t',
    usecols=["Gene","CNCCs_LFC_human_ref", "CNCCs_TPM_total", "CNCCs_gene_ase_type",
                "iPSCs_LFC_human_ref", "iPSCs_TPM_total", "iPSCs_gene_ase_type",
                "Pooled_mesenchyme_LFC_human_ref", "Pooled_mesenchyme_TPM_total", "Pooled_mesenchyme_gene_ase_type",
                "Pooled_prechondral_mesenchyme_LFC_human_ref", "Pooled_prechondral_mesenchyme_TPM_total", "Pooled_prechondral_mesenchyme_gene_ase_type",
                "9_Chondrocytes_LFC_human_ref", "9_Chondrocytes_TPM_total", "9_Chondrocytes_gene_ase_type"
             ])
merged_df = pd.merge(result_filtered, ase, left_on = 'Gene_symbol',right_on='Gene', how="left").drop("Gene", axis=1)

print(merged_df.head().to_string())

# merge with derived traits

#derived_traits = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/annoGenes_derTraits_chimp_skeletal.txt', header=0, sep='\t', usecols=["Gene", "Direction", "Derived traits - chimp", "Derived trait fraction - chimp", 'Direction trait matches (out of derived) - chimp', 'Fraction direction trait matches - chimp'])
derived_traits = pd.read_csv(f'/home/labs/davidgo/Collaboration/USEFUL DATASETS/Phenotypes/Neanderthal Chimp Phenotype annotation/annoGenes_derivedTraits/annoGenes_derTraits_chimp_skeleton.txt', header=0, sep='\t')

merged_df = pd.merge(merged_df, derived_traits, left_on = 'Gene_symbol',right_on='Gene', how="left").drop("Gene", axis=1)

print(merged_df.head().to_string())

# merge with KEgg
#kegg = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/KEGG_pathway_NEW.txt', header=0, sep='\t', usecols=["gene symbol","pathway name"])
kegg = pd.read_csv(f'/home/labs/davidgo/Collaboration/USEFUL DATASETS/Annotation DBs/KEGG/current/KEGG_pathway_NEW.txt', header=0, sep='\t', usecols=["gene symbol","pathway name"])
kegg_grouped = kegg.groupby('gene symbol')['pathway name'].agg(list).reset_index()

merged_df = pd.merge(merged_df, kegg_grouped, left_on = 'Gene_symbol',right_on='gene symbol', how="left").drop("gene symbol", axis=1)

print(merged_df.head().to_string())

# merge with gene cards
galc = pd.read_csv(f'/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/GeneCards/GeneAlaCart_gene_annotation/GALC_gene_annotation.txt', header=0, sep='\t')
skeletal_terms = pd.read_csv(f'/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/GeneCards/GeneCards skeletal terms/25Jun2024/input/GeneCards_skeletal_search.txt', header=0, sep='\t')

merged_df = pd.merge(merged_df, galc, left_on = 'Gene_symbol',right_on='Symbol', how="left").drop("Symbol", axis=1)
merged_df = pd.merge(merged_df, skeletal_terms, left_on = 'Gene_symbol',right_on='Symbol', how="left").drop("Symbol", axis=1)

print(merged_df.head().to_string())

# merge with tau
tau = pd.read_csv(f'/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/HPA/HPA_TAU.txt', usecols=["Gene","TAU score - Tissue", "TAU score - Single Cell Type"], header=0, sep='\t')
merged_df = pd.merge(merged_df, tau, left_on = 'Gene_symbol',right_on='Gene', how="left").drop("Gene", 1)
print(merged_df.head().to_string())

# merge with gene organizer
#geneorganizer = pd.read_excel('/home/labs/davidgo/Collaboration/humanMPRA/files_from_dropbox/ORGANizer_World_ALL_v14.xlsx', sheet_name="Systems", header=0)
geneorganizer = pd.read_excel('/home/labs/davidgo/Collaboration/USEFUL DATASETS/Annotation DBs/Gene ORGANizer/HPO 2024-04-26 Release/ORGANizer_World_HPO_ALL_HPO 2024-04-26 Release.xlsx', sheet_name="Systems", header=0)

def get_body_systems(row):
    body_systems = []
    for col in geneorganizer.columns:
        if col != 'gene_symbol' and col != 'entrez_ID' and row[col] == 1:
            body_systems.append(col)
    return body_systems
geneorganizer['Body Systems'] = geneorganizer.apply(get_body_systems, axis=1)
geneorganizer_filtered = geneorganizer[['gene_symbol','Body Systems']]

merged_df = pd.merge(merged_df, geneorganizer_filtered, left_on = 'Gene_symbol',right_on='gene_symbol', how="left").drop("gene_symbol", axis=1)
print(merged_df.head().to_string())

# merge with human and chimp expression data (genevive housman)

# save merged df
merged_df.to_csv(f'./top_candidates/chondrocytes/annotations_per_gene.csv', header=True, index=False)
