# this script takes the mpranalyze comparative input for oligo seq_325030_chr14:90968484-90968753_SCREEN_a2_L3 and create two new dataframe.
# in each of the new dataframes half of the columns is filled with 0 so that hopefully mpranalyze comparative will run faster (running on all the data reaches runtime limit of new_long)
# the column number (e.g barcode) that is filled with 0 has to be matched between the replicates for the ancestral allele, and also between the replicates for the derived allele (ut it does not have to match between ancestral and derived)
# I decided to fill columns with 0 instead of just taking half of the columns as this means that we can use the same depth factors easily and the same annot file

import pandas as pd
import numpy as np

# read in input for mprnalyze comparative and choose row that contains this oligo
df_full_rna = pd.read_csv(f'./chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header =0, index_col=0)
df_rna = df_full_rna.filter(items=['seq_325030_chr14:90968484-90968753_SCREEN_a2_L3'], axis=0)
print(np.count_nonzero(df_rna))
print(len(df_rna.columns))

df_full_dna = pd.read_csv(f'./chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header =0, index_col=0)
df_dna = df_full_dna.filter(items=['seq_325030_chr14:90968484-90968753_SCREEN_a2_L3'], axis=0)
print(np.count_nonzero(df_dna))
print(len(df_dna.columns))

# Extract unique numbers from column names for both ancestral and derived
numbers_ancestral = set(int(col.split('_')[-1]) for col in df_rna.columns if df_rna[col].iloc[0] != 0 and '_ancestral_' in col)
numbers_derived = set(int(col.split('_')[-1]) for col in df_rna.columns if df_rna[col].iloc[0] != 0 and '_derived_' in col)

# Create a random permutation of unique numbers for both ancestral and derived
random_numbers_ancestral = np.random.permutation(list(numbers_ancestral))
print(len(random_numbers_ancestral))
random_numbers_derived = np.random.permutation(list(numbers_derived))
print(len(random_numbers_derived))

# Select half of the numbers for both ancestral and derived
half_numbers_ancestral = random_numbers_ancestral[:len(random_numbers_ancestral)//2]
print(len(half_numbers_ancestral))
half_numbers_derived = random_numbers_derived[:len(random_numbers_derived)//2]
print(len(half_numbers_derived))

# Create masks for selecting ancestral and derived columns with 0s
mask_ancestral_group1 = df_rna.columns.str.endswith(tuple(f"_ancestral_{number}" for number in half_numbers_ancestral))
mask_derived_group1 = df_rna.columns.str.endswith(tuple(f"_derived_{number}" for number in half_numbers_derived))

# Combine masks to create final masks for group1
mask_group1 = np.logical_or(mask_ancestral_group1, mask_derived_group1)

#create masnk for group2
mask_group2 = ~mask_group1

# Create new DataFrames with selected columns
df_dna_group1 = df_dna.copy()
df_dna_group2 = df_dna.copy()

# Fill selected columns with 0s
df_dna_group1.loc[:, mask_group2] = 0.0
df_dna_group2.loc[:, mask_group1] = 0.0

# apply exactly the same masks also to rna input 
df_rna_group1 = df_rna.copy()
df_rna_group2 = df_rna.copy()

df_rna_group1.loc[:, mask_group2] = 0.0
df_rna_group2.loc[:, mask_group1] = 0.0

for df in [df_rna_group1, df_rna_group2, df_dna_group1, df_dna_group2]:
    print(np.count_nonzero(df))
    print((df == 0).sum().sum())
    print(len(df. columns))

# save dfs
df_dna_group1.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group1_non_zero_method.txt', sep='\t', header=False, index=True, index_label="locus")
df_dna_group2.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group2_non_zero_method.txt', sep='\t', header=False, index=True, index_label="locus")
df_rna_group1.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group1_non_zero_method.txt', sep='\t', header=False, index=True, index_label="locus")
df_rna_group2.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group2_non_zero_method.txt', sep='\t', header=False, index=True, index_label="locus")
