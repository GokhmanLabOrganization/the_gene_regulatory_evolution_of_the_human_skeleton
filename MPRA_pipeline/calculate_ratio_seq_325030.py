import pandas as pd
import numpy as np
import sys

# read in dataframe to take header from:
df_for_header = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header =0)
cols = df_for_header.columns.tolist()

# read in other dfs
dna_group1 = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group1.txt', sep='\t', header=None, names=cols)
dna_group2 = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group2.txt', sep='\t', header=None, names=cols)
rna_group1 = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group1.txt', sep='\t', header=None, names=cols)
rna_group2 = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group2.txt', sep='\t', header=None, names=cols)

for df in [dna_group1, dna_group2, rna_group1, rna_group2]:
    print((df == 0).sum().sum())

# # sum the ancestral and derived count for each df
# dna_group1_ancestral = dna_group1.loc[:, dna_group1.columns.str.contains('ancestral')].sum(axis=1).values[0]
# dna_group1_derived = dna_group1.loc[:, dna_group1.columns.str.contains('derived')].sum(axis=1).values[0]

# dna_group2_ancestral = dna_group2.loc[:, dna_group2.columns.str.contains('ancestral')].sum(axis=1).values[0]
# dna_group2_derived = dna_group2.loc[:, dna_group2.columns.str.contains('derived')].sum(axis=1).values[0]

# rna_group1_ancestral = rna_group1.loc[:, rna_group1.columns.str.contains('ancestral')].sum(axis=1).values[0]
# rna_group1_derived = rna_group1.loc[:, rna_group1.columns.str.contains('derived')].sum(axis=1).values[0]

# rna_group2_ancestral = rna_group2.loc[:, rna_group2.columns.str.contains('ancestral')].sum(axis=1).values[0]
# rna_group2_derived = rna_group2.loc[:, rna_group2.columns.str.contains('derived')].sum(axis=1).values[0]

# rna_list = [rna_group1_ancestral, rna_group1_derived, rna_group2_ancestral, rna_group2_derived]
# dna_list = [dna_group1_ancestral, dna_group1_derived, dna_group2_ancestral, dna_group2_derived]

# #chondro L3a2 - read counts for cpm normalization
# rep1DNA = 86012916
# rep1RNA = 262288284
# rep2DNA = 84852307
# rep2RNA = 267028259
# rep3DNA = 64816951
# rep3RNA = 252832900
# DNA_reads = rep1DNA+rep2DNA+rep3DNA
# RNA_reads = rep1RNA+rep2RNA+rep3RNA


# for i, rna in enumerate(rna_list):
    # rna_pseudo = rna + 1
    # dna_pseudo = dna[i] + 1
    # rna_cpm = (rna_pseudo*1000000)/RNA_reads
    # dna_cpm = (dna_pseudo*1000000)/DNA_reads
    # ratio = rna_cpm/dna_cpm
    # ratio_log = np.log2(ratio)
    # print(rna, rna_pseudo, dna, dna_pseudo, rna_cpm, dna_cpm, ratio, ratio_log)