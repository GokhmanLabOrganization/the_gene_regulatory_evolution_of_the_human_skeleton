import pandas as pd
import numpy as np

def count_nonzero_columns(group):
    return (np.count_nonzero(df))
    
df_with_names = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header =0)
cols = df_with_names.columns.tolist()

group1_dna = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group1.txt', sep='\t', header=None, names=cols)
group2_dna = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group2.txt', sep='\t', header=None, names=cols)
group1_rna = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group1.txt', sep='\t', header=None, names=cols)
group2_rna = pd.read_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/chondrocytes/L3a2/output/mpranalyze_comparative/chondrocytes_L3a2_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_seq_325030_group2.txt', sep='\t', header=None, names=cols)

for df in [group1_dna, group2_dna, group1_rna, group2_rna]:

    # Initialize counters
    rep1_count_ancestral = count_nonzero_columns(df.filter(like='rep1_ancestral'))
    rep2_count_ancestral = count_nonzero_columns(df.filter(like='rep2_ancestral'))
    rep3_count_ancestral = count_nonzero_columns(df.filter(like='rep3_ancestral'))
    rep1_count_derived = count_nonzero_columns(df.filter(like='rep1_derived'))
    rep2_count_derived = count_nonzero_columns(df.filter(like='rep2_derived'))
    rep3_count_derived = count_nonzero_columns(df.filter(like='rep3_derived'))


    # Print the results
    print("rep1 count ancestral:", rep1_count_ancestral)
    print("rep2 count ancestral:", rep2_count_ancestral)
    print("rep3 count ancestral:", rep3_count_ancestral)
    print("rep1 count derived:", rep1_count_derived)
    print("rep2 count derived:", rep2_count_derived)
    print("rep3 count derived:", rep3_count_derived)