import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from ast import literal_eval


# read in comb df for L1a2
L1a2 = pd.read_csv('./chondrocytes/L1a2/output/activity_after_filter/comb_df_adjusted_fdr.csv', converters={"RNA_chondrocytes_rep1": literal_eval, "RNA_chondrocytes_rep2": literal_eval, "RNA_chondrocytes_rep3": literal_eval, "DNA_chondrocytes_rep1": literal_eval, "DNA_chondrocytes_rep2": literal_eval, "DNA_chondrocytes_rep3": literal_eval})

# filter for the relevant sequences
L1a2_filtered = L1a2[L1a2['oligo'].str.contains('seq_59045_chr12:108511253-108511522_SCREEN_ancestral_a2_L1|seq_79469_chr3:197439120-197439389_SCREEN_ancestral_a2_L1|seq_84596_chr3:14692942-14693211_SCREEN_derived_a2_L1|seq_64068_chr7:50599046-50599315_SCREEN_ancestral_a2_L1|scrambled')].copy()
# L1a2_filtered.to_csv('./additional/saturation_oligos/L1a2_filtered.csv', header=True, index = False)


# sum up the UMI count lists of the different replicates for RNA
sum_lists_RNA = lambda row: [sum(x) for x in zip(row['RNA_chondrocytes_rep1'], row['RNA_chondrocytes_rep2'], row['RNA_chondrocytes_rep3'])]
L1a2_filtered['sum_list_RNA'] = L1a2_filtered.apply(sum_lists_RNA, axis=1)
# sum up the UMI count lists of the different replicates for DNA
sum_lists_DNA = lambda row: [sum(x) for x in zip(row['DNA_chondrocytes_rep1'], row['DNA_chondrocytes_rep2'], row['DNA_chondrocytes_rep3'])]
L1a2_filtered['sum_list_DNA'] = L1a2_filtered.apply(sum_lists_DNA, axis=1)

# add pseudocounts to DNA
add_one_to_list = lambda x: [i + 1 for i in x]
L1a2_filtered['sum_list_DNA_pseudocounts'] = L1a2_filtered['sum_list_DNA'].apply(add_one_to_list)


# divide RNA by DNA
RNA_DNA_ratio_per_barcode = lambda row: [a / b for a, b in zip(row['sum_list_RNA'], row['sum_list_DNA_pseudocounts'])]
L1a2_filtered['RNA_DNA_ratio_list'] = L1a2_filtered.apply(RNA_DNA_ratio_per_barcode, axis=1)

# scrambled
L1a2_scrambled = L1a2_filtered[L1a2_filtered['oligo'].str.contains('scrambled')].copy()
combined_list_scrambled = L1a2_scrambled['RNA_DNA_ratio_list'].sum()
scrambled_combined_df = pd.DataFrame({'RNA_DNA_ratio_list': [combined_list_scrambled], 'oligo': ['scrambled']})

L1a2_filtered_wo_scrambled = L1a2_filtered[~L1a2_filtered['oligo'].str.contains('scrambled')]
L1a2_with_scrambled = pd.concat([L1a2_filtered_wo_scrambled, scrambled_combined_df], join="inner")


# save df for qc
L1a2_with_scrambled.to_csv('./additional/saturation_oligos/L1a2_ratio_per_barcode_with_scrambled.csv', header=True, index = False)

# plot a histogramm 
plt.clf()
fig, axs = plt.subplots(nrows=len(L1a2_with_scrambled), ncols=1, sharex=True)
for i, (index, row) in enumerate(L1a2_with_scrambled.iterrows()):
    axs[i].hist(row['RNA_DNA_ratio_list'], bins=60, range=(0,120), label=f'{row["oligo"]}')
    axs[i].legend()

plt.xlabel('RNA/DNA ratio')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig(f'./additional/saturation_oligos/hist_rna_dna_barcode_modified.pdf')
plt.savefig(f'./additional/saturation_oligos/hist_rna_dna_barcode_modified.png', dpi=330)

plt.clf()
fig, axs = plt.subplots(nrows=len(L1a2_with_scrambled), ncols=1, sharex=True)
for i, (index, row) in enumerate(L1a2_with_scrambled.iterrows()):
    axs[i].hist(row['RNA_DNA_ratio_list'], label=f'{row["oligo"]}')
    axs[i].legend()

plt.xlabel('RNA/DNA ratio')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig(f'./additional/saturation_oligos/hist_rna_dna_barcode.pdf')
plt.savefig(f'./additional/saturation_oligos/hist_rna_dna_barcode.png', dpi=330)