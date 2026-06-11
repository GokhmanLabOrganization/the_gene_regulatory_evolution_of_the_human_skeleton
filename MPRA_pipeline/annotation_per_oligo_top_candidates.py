# this script merges summary table fomr mpra, output from yardens script (Summary_PerRegion2AllGenes.tsv) and output from EAT to find top candidates
import pandas as pd
import sys
import pyBigWig

cells = sys.argv[1]

mpra = pd.read_csv(f'./{cells}/comparative_analysis_combined/humanMPRA_with_seq_final2.csv', header=0, usecols=["oligo","RNA_DNA_ratio_log2_ancestral", "RNA_DNA_ratio_log2_derived", "DNA_counts_raw_ancestral","DNA_counts_raw_derived","barcode_count_ancestral", "barcode_count_derived", "normalized_activity_estimate_ancestral", "normalized_activity_estimate_derived", "activity_fdr_ancestral", "activity_fdr_derived", "activity_ancestral", "activity_derived", "logFC_derived_vs_ancestral", "differential_activity_fdr", "differential_activity"])
genes = pd.read_csv(f'./region_gene_link/{cells}/yardens_script/output_Nadav2_16-12-2024/Summary_PerRegion2AllGenes.tsv', sep='\t', header=0) #NM 18-12-2024 Changed the output to updated output including L4a1
eat = pd.read_excel(f'./EAT/Nadav_final_2024-12-20_18-10-47.xlsx', header=0) #NM 18-12-2024 Updated to include a new version with L4a1
variant_count = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA_raw/sequences.xlsx', sheet_name = "oligos", header=0)
variant_count_L4 = pd.read_excel(f'/home/labs/davidgo/Collaboration/humanMPRA_raw/L4/L4_v5.xlsx', sheet_name = "test_oligos", header=0) #Added the L4 library
variant_count_L4_a1 = variant_count_L4[variant_count_L4['adapter'] == 'a1']


# create id to merge on
mpra['id'] = mpra['oligo'].str.split('_').str[2]

genes['id'] = genes['chromosome'].astype(str) + ':' + (genes['start'] + 1).astype(str) + '-' + genes['end'].astype(str)

eat['id'] = eat['chr'].astype(str) + ':' + (eat['start'] + 1).astype(str) + '-' + eat['end'].astype(str)
eat.drop(columns=['chr', 'start', 'end'], inplace=True)

variant_count['id'] = variant_count['chr'].astype(str) + ':' + variant_count['start'].astype(str) + '-' + variant_count['end'].astype(str)

variant_count_L4_a1['id'] = variant_count_L4_a1['chr'].astype(str) + ':' + variant_count_L4_a1['start'].astype(str) + '-' + variant_count_L4_a1['end'].astype(str)

# merge variant count
variant_count_filtered = variant_count[["id", "variants_count", "variants_genomic"]]
variant_count_L4_a1_filtered = variant_count_L4_a1[["id", "variants_count", "variants_genomic"]]
combined_variant_count_filtered = pd.concat([variant_count_filtered, variant_count_L4_a1_filtered], ignore_index=True)

merged_df = pd.merge(mpra, combined_variant_count_filtered, on='id', how="left")

# merge genes
merged_df = pd.merge(merged_df, genes, on='id', how="left")


# merge eat
merged_df = pd.merge(merged_df, eat, on='id', how="left")

# merge SCREEN annotation

# # merge phylop
# bw=pyBigWig.open("/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/phyloP_hg19/hg19.100way.phyloP100way.bw")
# def extract_values(row):
    # chr_val = str(row['chromosome'])  # Convert chromosome to string
    # start_val = int(row['start'])
    # end_val = int(row['end'])
    
    # # Calculate average value
    # average_value = bw.stats(chr_val, start_val, end_val)[0]
    
    # # Calculate maximum value
    # max_value = bw.stats(chr_val, start_val, end_val, type="max")[0]
    
    # return pd.Series({'average_value': average_value, 'max_value': max_value})

# merged_df[['phylop_average', 'phylop_max']] = merged_df.apply(extract_values, axis=1)

# bw.close()

print(len(merged_df))
df_cleaned = merged_df.drop_duplicates()
print(len(df_cleaned))

# save df
df_cleaned.to_csv(f'./top_candidates/chondrocytes/humanMPRA_annotations_v3.csv', header=True, index = False)
