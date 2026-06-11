# this script performs some analysis that compares then neurons to the chondrocytes to better understand the issues with the neurons
# 1st analysis: comparing activity levels of negative controls between chondorcuyes and neurons
# 2nd analysis: compare DNA and RNA counts per barcode
# 3rd analysis: comparing distributions of rna/dna, rna (normalized to cpm and number of barcodes) and dna (normalized to cpm and number of barcodes)
# 4th analysis: compare number of barcodes per oligo
# 5th analysis: correlations of counts per barcode between chondrocytes and neurons (rna-rna, dna-dna, dna-rna, rna-dna)
# analysis 6: correlations of counts per barcode and ratio per barcode between replicates (neurons vs chondrocytes)
# analysis 7: correlations between dna and rna counts per barcode (for chondrocytes and neurons) - to compare to figure 5
# analysis 8: histogram of rna counts per oligo comparing neurons and chondrocytes + scatter (counts normalized to library size)
# analysis 9: scatter comparing ratio per oligo in chondrocytes vs neurons
# analysis 10: scatter comparing ratio per oligo in chondrocytes vs neurons - only taking oligos with DNA > 500

import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from scipy import stats


library = sys.argv[1]
adaptor = sys.argv[2]

comb_df_path_neurons = f'./neurons/{library}{adaptor}/output/activity_after_filter/comb_df_adjusted_fdr.csv'
comb_df_path_chondro = f'./chondrocytes/{library}{adaptor}/output/activity_after_filter/comb_df_adjusted_fdr.csv'

neuron_df = pd.read_csv(comb_df_path_neurons)
chondro_df = pd.read_csv(comb_df_path_chondro)


# # 1st analysis: comparing activity levels of negative controls between chondorcuyes and neurons

# filtered_neuron_df = neuron_df[neuron_df['control_type'].isin(['NegCtrl_not_active', 'scrambled_control', 'NegCtrl_non_SCREEN'])]
# filtered_chondro_df = chondro_df[chondro_df['control_type'].isin(['NegCtrl_not_active', 'scrambled_control', 'NegCtrl_non_SCREEN'])]

# # boxplot
# neuro_chondro_concat_df = pd.concat([filtered_neuron_df.assign(dataset='neurons'), filtered_chondro_df.assign(dataset='chondrocytes')])

# plt.clf()
# plt.figure()
# sns.boxplot(x='control_type', y='ratio_log_rep_comb', hue='dataset', data=neuro_chondro_concat_df)
# plt.savefig(f'./additional/chondro_neuron_comparison/ratio_negative_controls_boxplot_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/ratio_negative_controls_boxplot_{library}{adaptor}.png', dpi=330)

# # scatter
# neuro_chondro_merge_df = filtered_neuron_df.merge(filtered_chondro_df, on="oligo", suffixes=('_neuron', '_chondro'))

# plt.clf()
# sns.scatterplot(data = neuro_chondro_merge_df, x = f"ratio_log_rep_comb_neuron", y = f"ratio_log_rep_comb_chondro")
# plt.axes().set_aspect('equal')
# plt.savefig(f'./additional/chondro_neuron_comparison/ratio_negative_controls_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/ratio_negative_controls_scatter_{library}{adaptor}.png', dpi=330)

# # 2nd analysis: compare DNA and RNA counts per barcode

umi_chondro = pd.read_csv(f'./chondrocytes/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t', header=0)
umi_neuron = pd.read_csv(f'./neurons/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t', header=0)

# umi_chondro['RNA'] = umi_chondro[[col for col in umi_chondro.columns if 'RNA' in col]].sum(axis=1)
# umi_chondro['DNA'] = umi_chondro[[col for col in umi_chondro.columns if 'DNA' in col]].sum(axis=1)

# umi_neuron['RNA'] = umi_neuron[[col for col in umi_neuron.columns if 'RNA' in col]].sum(axis=1)
# umi_neuron['DNA'] = umi_neuron[[col for col in umi_neuron.columns if 'DNA' in col]].sum(axis=1)

# umi_chondro['rna_normalized'] = (umi_chondro['RNA'] / umi_chondro['RNA'].sum()) * 1e6
# umi_chondro['dna_normalized'] = (umi_chondro['DNA'] / umi_chondro['DNA'].sum()) * 1e6

# umi_neuron['rna_normalized'] = (umi_neuron['RNA'] / umi_neuron['RNA'].sum()) * 1e6
# umi_neuron['dna_normalized'] = (umi_neuron['DNA'] / umi_neuron['DNA'].sum()) * 1e6

# plt.clf()
# fig, axs = plt.subplots(2, 2, sharex="col")

# axs[0, 0].hist(umi_neuron['rna_normalized'], bins=100, range=[0, 5])
# axs[0, 1].hist(umi_neuron['dna_normalized'], bins=100, range=[0, 1])
# axs[1, 0].hist(umi_chondro['rna_normalized'], bins=100, range=[0, 5])
# axs[1, 1].hist(umi_chondro['dna_normalized'], bins=100, range=[0, 1])

# axs[0, 0].set_title('Neurons - RNA Histogram')
# axs[0, 1].set_title('Neurons - DNA Histogram')
# axs[1, 0].set_title('Chondrocytes - RNA Histogram')
# axs[1, 1].set_title('Chondrocytes - DNA Histogram')

# plt.savefig(f'./additional/chondro_neuron_comparison/counts_per_barcode_hist_xlim_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/counts_per_barcode_hist_xlim_{library}{adaptor}.png', dpi=330)

# 3rd analysis: comparing distributions of rna/dna, rna (normalized to cpm and number of barcodes) and dna (normalized to cpm and number of barcodes)

# cpm_neurons = pd.read_csv(f'./neurons/{library}{adaptor}/output/activity_after_filter/ratio_after_filter.csv')
# cpm_chondro = pd.read_csv(f'./chondrocytes/{library}{adaptor}/output/activity_after_filter/ratio_after_filter.csv')

# cpm_neurons["RNA_cpm_bc_norm"] = cpm_neurons["RNA_cpm_rep_comb"]/cpm_neurons["count_rep_comb"]
# cpm_neurons["DNA_cpm_bc_norm"] = cpm_neurons["DNA_cpm_rep_comb"]/cpm_neurons["count_rep_comb"]
# cpm_chondro["RNA_cpm_bc_norm"] = cpm_chondro["RNA_cpm_rep_comb"]/cpm_chondro["count_rep_comb"]
# cpm_chondro["DNA_cpm_bc_norm"] = cpm_chondro["DNA_cpm_rep_comb"]/cpm_chondro["count_rep_comb"]

# plt.clf()
# fig, axs = plt.subplots(2, 3, sharex="col", figsize=[9, 4.8])
# sns.histplot(ax=axs[0, 0], data=neuron_df, x=f"ratio_log_rep_comb", hue = "activity", hue_order = ["non_active", "active"], multiple="stack", palette=["darkgrey", "darkorange"], edgecolor='none', alpha=1,binwidth=0.03 )
# sns.histplot(ax=axs[1, 0], data=chondro_df, x=f"ratio_log_rep_comb", hue = "activity", hue_order = ["non_active", "active"], multiple="stack", palette=["darkgrey", "darkorange"], edgecolor='none', alpha=1,binwidth=0.03 )

# sns.histplot(ax=axs[0, 1], data=cpm_neurons, x=f"RNA_cpm_bc_norm", edgecolor='none')
# sns.histplot(ax=axs[1, 1], data=cpm_chondro, x=f"RNA_cpm_bc_norm", edgecolor='none')

# sns.histplot(ax=axs[0, 2], data=cpm_neurons, x=f"DNA_cpm_bc_norm", edgecolor='none')
# sns.histplot(ax=axs[1, 2], data=cpm_chondro, x=f"DNA_cpm_bc_norm", edgecolor='none')

# fig.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_dna_dist_bc_normalized_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_dna_dist_bc_normalized_{library}{adaptor}.png', dpi=330)

# plt.clf()
# fig, axs = plt.subplots(2, 3, sharex="col", figsize=[9, 4.8])
# sns.histplot(ax=axs[0, 0], data=neuron_df, x=f"ratio_log_rep_comb", hue = "activity", hue_order = ["non_active", "active"], multiple="stack", palette=["darkgrey", "darkorange"], edgecolor='none', alpha=1,binwidth=0.03 )
# sns.histplot(ax=axs[1, 0], data=chondro_df, x=f"ratio_log_rep_comb", hue = "activity", hue_order = ["non_active", "active"], multiple="stack", palette=["darkgrey", "darkorange"], edgecolor='none', alpha=1,binwidth=0.03 )

# sns.histplot(ax=axs[0, 1], data=cpm_neurons, x=f"RNA_cpm_bc_norm", edgecolor='none')
# sns.histplot(ax=axs[1, 1], data=cpm_chondro, x=f"RNA_cpm_bc_norm", edgecolor='none')
# axs[0, 1].set_xlim([0, 0.2])
# axs[1, 1].set_xlim([0, 0.2])

# sns.histplot(ax=axs[0, 2], data=cpm_neurons, x=f"DNA_cpm_bc_norm", edgecolor='none')
# sns.histplot(ax=axs[1, 2], data=cpm_chondro, x=f"DNA_cpm_bc_norm", edgecolor='none')
# axs[0, 2].set_xlim([0, 0.15])
# axs[1, 2].set_xlim([0, 0.15])

# fig.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_dna_dist_bc_normalized_xlim_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_dna_dist_bc_normalized_xlim_{library}{adaptor}.png', dpi=330)

# # 4th analysis: compare number of barcodes per oligo

# plt.clf()
# fig, axs = plt.subplots(2, 6, sharex="col", figsize=[12, 4.8])
# sns.histplot(ax=axs[0, 0], data=neuron_df, x=f"bcs_RNA_rep1", edgecolor='none')
# sns.histplot(ax=axs[0, 1], data=neuron_df, x=f"bcs_RNA_rep2", edgecolor='none')
# sns.histplot(ax=axs[0, 2], data=neuron_df, x=f"bcs_RNA_rep3", edgecolor='none')
# sns.histplot(ax=axs[0, 3], data=neuron_df, x=f"bcs_DNA_rep1", edgecolor='none')
# sns.histplot(ax=axs[0, 4], data=neuron_df, x=f"bcs_DNA_rep2", edgecolor='none')
# sns.histplot(ax=axs[0, 5], data=neuron_df, x=f"bcs_DNA_rep3", edgecolor='none')
# sns.histplot(ax=axs[1, 0], data=chondro_df, x=f"bcs_RNA_rep1", edgecolor='none')
# sns.histplot(ax=axs[1, 1], data=chondro_df, x=f"bcs_RNA_rep2", edgecolor='none')
# sns.histplot(ax=axs[1, 2], data=chondro_df, x=f"bcs_RNA_rep3", edgecolor='none')
# sns.histplot(ax=axs[1, 3], data=chondro_df, x=f"bcs_DNA_rep1", edgecolor='none')
# sns.histplot(ax=axs[1, 4], data=chondro_df, x=f"bcs_DNA_rep2", edgecolor='none')
# sns.histplot(ax=axs[1, 5], data=chondro_df, x=f"bcs_DNA_rep3", edgecolor='none')
# fig.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_per_oligo_distribution_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_per_oligo_distribution_{library}{adaptor}.png', dpi=330)

# # scatter
merged_all = neuron_df.merge(chondro_df, on="oligo", suffixes=('_neuron', '_chondro'))
# plt.clf()
# fig, axs = plt.subplots(2, 3, sharex="row", sharey="row", figsize=[9, 4.8])
# sns.scatterplot(ax=axs[0, 0], data=merged_all, x=f"bcs_RNA_rep1_neuron", y="bcs_RNA_rep1_chondro", s=1, linewidth=0)
# sns.scatterplot(ax=axs[0, 1], data=merged_all, x=f"bcs_RNA_rep2_neuron", y="bcs_RNA_rep2_chondro", s=1, linewidth=0)
# sns.scatterplot(ax=axs[0, 2], data=merged_all, x=f"bcs_RNA_rep3_neuron", y="bcs_RNA_rep3_chondro", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 0], data=merged_all, x=f"bcs_DNA_rep1_neuron", y="bcs_DNA_rep1_chondro", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 1], data=merged_all, x=f"bcs_DNA_rep2_neuron", y="bcs_DNA_rep2_chondro", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 2], data=merged_all, x=f"bcs_DNA_rep3_neuron", y="bcs_DNA_rep3_chondro", s=1, linewidth=0)
# axs[0, 0].set_aspect('equal')
# axs[0, 1].set_aspect('equal')
# axs[0, 2].set_aspect('equal')
# axs[1, 0].set_aspect('equal')
# axs[1, 1].set_aspect('equal')
# axs[1, 2].set_aspect('equal')
# axs[0, 0].yaxis.get_label().set_visible(True)
# axs[0, 1].yaxis.get_label().set_visible(True)
# axs[0, 2].yaxis.get_label().set_visible(True)
# axs[1, 0].yaxis.get_label().set_visible(True)
# axs[1, 1].yaxis.get_label().set_visible(True)
# axs[1, 2].yaxis.get_label().set_visible(True)
# fig.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_per_oligo_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_per_oligo_scatter_{library}{adaptor}.png', dpi=330)

# 5th analysis: correlations of counts per barcode between chondrocytes and neurons (rna-rna, dna-dna, dna-rna, rna-dna)
# print(umi_chondro['RNA_chondrocytes_rep1'].sum())

umi_chondro['rna_rep1_normalized'] = (umi_chondro['RNA_chondrocytes_rep1'] / umi_chondro['RNA_chondrocytes_rep1'].sum()) * 1e6
umi_chondro['dna_rep1_normalized'] = (umi_chondro['DNA_chondrocytes_rep1'] / umi_chondro['DNA_chondrocytes_rep1'].sum()) * 1e6
umi_chondro['rna_rep2_normalized'] = (umi_chondro['RNA_chondrocytes_rep2'] / umi_chondro['RNA_chondrocytes_rep2'].sum()) * 1e6
umi_chondro['dna_rep2_normalized'] = (umi_chondro['DNA_chondrocytes_rep2'] / umi_chondro['DNA_chondrocytes_rep2'].sum()) * 1e6
umi_chondro['rna_rep3_normalized'] = (umi_chondro['RNA_chondrocytes_rep3'] / umi_chondro['RNA_chondrocytes_rep3'].sum()) * 1e6
umi_chondro['dna_rep3_normalized'] = (umi_chondro['DNA_chondrocytes_rep3'] / umi_chondro['DNA_chondrocytes_rep3'].sum()) * 1e6

umi_neuron['rna_rep1_normalized'] = (umi_neuron['RNA_neurons_rep1'] / umi_neuron['RNA_neurons_rep1'].sum()) * 1e6
umi_neuron['dna_rep1_normalized'] = (umi_neuron['DNA_neurons_rep1'] / umi_neuron['DNA_neurons_rep1'].sum()) * 1e6
umi_neuron['rna_rep2_normalized'] = (umi_neuron['RNA_neurons_rep2'] / umi_neuron['RNA_neurons_rep2'].sum()) * 1e6
umi_neuron['dna_rep2_normalized'] = (umi_neuron['DNA_neurons_rep2'] / umi_neuron['DNA_neurons_rep2'].sum()) * 1e6
umi_neuron['rna_rep3_normalized'] = (umi_neuron['RNA_neurons_rep3'] / umi_neuron['RNA_neurons_rep3'].sum()) * 1e6
umi_neuron['dna_rep3_normalized'] = (umi_neuron['DNA_neurons_rep3'] / umi_neuron['DNA_neurons_rep3'].sum()) * 1e6

# umi_merged = umi_neuron.merge(umi_chondro, on="oligo_bc", suffixes=('_neuron', '_chondro'))
# plt.clf()
# fig, axs = plt.subplots(2, 2)
# sns.scatterplot(ax=axs[0, 0], data=umi_merged, x=f"dna_rep1_normalized_chondro", y="dna_rep1_normalized_neuron", s=1, linewidth=0)
# sns.scatterplot(ax=axs[0, 1], data=umi_merged, x=f"rna_rep1_normalized_chondro", y="dna_rep1_normalized_neuron", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 0], data=umi_merged, x=f"dna_rep1_normalized_chondro", y="rna_rep1_normalized_neuron", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 1], data=umi_merged, x=f"rna_rep1_normalized_chondro", y="rna_rep1_normalized_neuron", s=1, linewidth=0)
# fig.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_counts_neurons_vs_chondro_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_counts_neurons_vs_chondro_scatter_{library}{adaptor}.png', dpi=330)

# d_c_d_n = umi_merged['dna_rep1_normalized_chondro'].corr(umi_merged["dna_rep1_normalized_neuron"])
# r_c_d_n = umi_merged['rna_rep1_normalized_chondro'].corr(umi_merged["dna_rep1_normalized_neuron"])
# d_c_r_n = umi_merged['dna_rep1_normalized_chondro'].corr(umi_merged["rna_rep1_normalized_neuron"])
# r_c_r_n = umi_merged['rna_rep1_normalized_chondro'].corr(umi_merged["rna_rep1_normalized_neuron"])

# print(f'corr dna(chondrocytes)- dna(neurons) {d_c_d_n}')
# print(f'corr rna(chondrocytes)- dna(neurons) {r_c_d_n}')
# print(f'corr dna(chondrocytes)- rna(neurons) {d_c_r_n}')
# print(f'corr rna(chondrocytes)- rna(neurons) {r_c_r_n}')

# # analysis 6: correlations of counts per barcode and ratio per barcode between replicates (neurons vs chondrocytes)

# umi_neuron["ratio_rep1"] = np.log2(umi_neuron["rna_rep1_normalized"]/umi_neuron["dna_rep1_normalized"])
# umi_neuron["ratio_rep2"] = np.log2(umi_neuron["rna_rep2_normalized"]/umi_neuron["dna_rep2_normalized"])

# umi_chondro["ratio_rep1"] = np.log2(umi_chondro["rna_rep1_normalized"]/umi_chondro["dna_rep1_normalized"])
# umi_chondro["ratio_rep2"] = np.log2(umi_chondro["rna_rep2_normalized"]/umi_chondro["dna_rep2_normalized"])

# plt.clf()
# fig, axs = plt.subplots(2, 3)
# sns.scatterplot(ax=axs[0, 0], data=umi_neuron, x=f"dna_rep1_normalized", y="dna_rep2_normalized", s=1, linewidth=0)
# sns.scatterplot(ax=axs[0, 1], data=umi_neuron, x=f"rna_rep1_normalized", y="rna_rep2_normalized", s=1, linewidth=0)
# sns.scatterplot(ax=axs[0, 2], data=umi_neuron, x=f"ratio_rep1", y="ratio_rep2", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 0], data=umi_chondro, x=f"dna_rep1_normalized", y="dna_rep2_normalized", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 1], data=umi_chondro, x=f"rna_rep1_normalized", y="rna_rep2_normalized", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1, 2], data=umi_chondro, x=f"ratio_rep1", y="ratio_rep2", s=1, linewidth=0)
# axs[0, 0].set_aspect('equal')
# axs[0, 1].set_aspect('equal')
# axs[0, 2].set_aspect('equal')
# axs[1, 0].set_aspect('equal')
# axs[1, 1].set_aspect('equal')
# axs[1, 2].set_aspect('equal')
# fig.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_counts_corr_between_reps_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/barcode_counts_corr_between_reps_scatter_{library}{adaptor}.png', dpi=330)

# print(umi_chondro["ratio_rep2"].head())
# print(umi_chondro["ratio_rep1"].head())

# umi_chondro_no_inf = umi_chondro.replace([np.inf, -np.inf], np.nan)
# umi_neuron_no_inf = umi_neuron.replace([np.inf, -np.inf], np.nan)

# umi_chondro_no_inf.dropna(subset=["ratio_rep1", "ratio_rep2"], inplace=True)
# umi_neuron_no_inf.dropna(subset=["ratio_rep1", "ratio_rep2"], inplace=True)

# c_dna = umi_chondro['dna_rep1_normalized'].corr(umi_chondro["dna_rep2_normalized"])
# c_rna = umi_chondro['rna_rep1_normalized'].corr(umi_chondro["rna_rep2_normalized"])
# c_ratio = umi_chondro_no_inf['ratio_rep1'].corr(umi_chondro_no_inf["ratio_rep2"])

# n_dna = umi_neuron['dna_rep1_normalized'].corr(umi_neuron["dna_rep2_normalized"])
# n_rna = umi_neuron['rna_rep1_normalized'].corr(umi_neuron["rna_rep2_normalized"])
# n_ratio = umi_neuron_no_inf['ratio_rep1'].corr(umi_neuron_no_inf["ratio_rep2"])

# print("chondrocytes", c_dna, c_rna, c_ratio)
# print("neurons", n_dna, n_rna, n_ratio)

# # analysis 7: correlations between dna and rna counts per barcode (for chondrocytes and neurons) - to compare to figure 5

# plt.clf()
# fig, axs = plt.subplots(1, 2)
# values_neurons = np.vstack([umi_neuron[f"dna_rep1_normalized"], umi_neuron[f"rna_rep1_normalized"]])
# kernel_neurons = stats.gaussian_kde(values_neurons)(values_neurons)
# values_chondro = np.vstack([umi_chondro[f"dna_rep1_normalized"], umi_chondro[f"rna_rep1_normalized"]])
# kernel_chondro = stats.gaussian_kde(values_chondro)(values_chondro)
# sns.scatterplot(ax=axs[0], data=umi_neuron, x=f"dna_rep1_normalized", y="rna_rep1_normalized", s=1, linewidth=0, c=kernel_neurons)
# sns.scatterplot(ax=axs[1], data=umi_chondro, x=f"dna_rep1_normalized", y="rna_rep1_normalized", s=1, linewidth=0, c=kernel_chondro)
# fig.tight_layout()
# axs[0].set_title('Neurons')
# axs[1].set_title('Chondrocytes')
# axs[0].set_aspect('equal')
# axs[1].set_aspect('equal')
# plt.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_vs_dna_barcode_counts_scatter_{library}{adaptor}_color.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_vs_dna_barcode_counts_scatter_{library}{adaptor}_color.png', dpi=330)

# corr_neurons = umi_neuron['dna_rep1_normalized'].corr(umi_neuron["rna_rep1_normalized"])
# corr_chondro = umi_chondro['dna_rep1_normalized'].corr(umi_chondro["rna_rep1_normalized"])

# print(f'corr dna(chondrocytes)- rna(chondrocytes) {corr_chondro}')
# print(f'corr dna(neurons)- rna(neurons) {corr_neurons}')

# # analysis 8: histogram of rna counts per oligo comparing neurons and chondrocytes + scatter (counts normalized to library size)

# for rep in ["RNA_rep1","RNA_rep2","RNA_rep3","RNA_rep_comb"]:
    # neuron_df[f'{rep}_normalized'] = (neuron_df[f'{rep}'] / neuron_df[f'{rep}'].sum()) * 1e6
    # chondro_df[f'{rep}_normalized'] = (chondro_df[f'{rep}'] / chondro_df[f'{rep}'].sum()) * 1e6


# plt.clf()
# fig, axes = plt.subplots(2,3, sharey=True, sharex=True)
# for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    # curr_df_neurons = neuron_df.copy()
    # curr_df_neurons[f"RNA_{rep}_normalized"] = curr_df_neurons[f"RNA_{rep}_normalized"].clip(upper=125)
    # curr_df_chondro = chondro_df.copy()
    # curr_df_chondro[f"RNA_{rep}_normalized"] = curr_df_chondro[f"RNA_{rep}_normalized"].clip(upper=125)
    # sns.histplot(ax=axes[0, n], data=curr_df_neurons, x=f"RNA_{rep}_normalized", edgecolor='none')
    # sns.histplot(ax=axes[1, n], data=curr_df_chondro, x=f"RNA_{rep}_normalized", edgecolor='none')
# plt.suptitle("Distribution RNA (cpm) - clipped to 125")
# plt.savefig(f'./additional/chondro_neuron_comparison/RNA_distribution_cpm_clipped125_{library}{adaptor}.png', transparent=True, dpi=330)
# plt.savefig(f'./additional/chondro_neuron_comparison/RNA_distribution_cpm_clipped125_{library}{adaptor}.pdf')

merged_all_with_normalized = neuron_df.merge(chondro_df, on="oligo", suffixes=('_neuron', '_chondro'))

# merged_all_with_normalized_no_na_rna = merged_all_with_normalized.dropna(subset=['RNA_rep_comb_normalized_neuron', 'RNA_rep_comb_normalized_chondro'], how='any')

# plt.clf()
# values = np.vstack([merged_all_with_normalized_no_na_rna[f"RNA_rep_comb_normalized_neuron"], merged_all_with_normalized_no_na_rna[f"RNA_rep_comb_normalized_chondro"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_all_with_normalized_no_na_rna, x = f"RNA_rep_comb_normalized_neuron", y = f"RNA_rep_comb_normalized_chondro", s=1, c=kernel)
# plt.axes().set_aspect('equal')
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_chondro_vs_neurons_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison/rna_chondro_vs_neurons_scatter_{library}{adaptor}.png', dpi=330)


# analysis 9: scatter comparing ratio per oligo in chondrocytes vs neurons

merged_all_with_normalized_no_na_ratio = merged_all_with_normalized.dropna(subset=['ratio_log_rep_comb_neuron', 'ratio_log_rep_comb_chondro'], how='any')

plt.clf()
values = np.vstack([merged_all_with_normalized_no_na_ratio[f"ratio_log_rep_comb_neuron"], merged_all_with_normalized_no_na_ratio[f"ratio_log_rep_comb_chondro"]])
kernel = stats.gaussian_kde(values)(values)
sns.scatterplot(data = merged_all_with_normalized_no_na_ratio, x = f"ratio_log_rep_comb_neuron", y = f"ratio_log_rep_comb_chondro", s=1, c=kernel)
plt.axes().set_aspect('equal')
plt.savefig(f'./additional/chondro_neuron_comparison/ratio_chondro_vs_neurons_scatter_{library}{adaptor}.pdf')
plt.savefig(f'./additional/chondro_neuron_comparison/ratio_chondro_vs_neurons_scatter_{library}{adaptor}.png', dpi=330)

# analysis 10: scatter comparing ratio per oligo in chondrocytes vs neurons - only taking oligos with DNA > 500
# recalculate normalization after filtering out oligos

merged_all_500 = merged_all[merged_all["DNA_rep_comb_neuron"]>500].copy()
merged_all_500[f'DNA_rep_comb_chondro_normalized500'] = (merged_all_500[f'DNA_rep_comb_chondro'] / merged_all_500[f'DNA_rep_comb_chondro'].sum()) * 1e6
merged_all_500[f'RNA_rep_comb_chondro_normalized500'] = (merged_all_500[f'RNA_rep_comb_chondro'] / merged_all_500[f'RNA_rep_comb_chondro'].sum()) * 1e6
merged_all_500[f'DNA_rep_comb_neuron_normalized500'] = (merged_all_500[f'DNA_rep_comb_neuron'] / merged_all_500[f'DNA_rep_comb_neuron'].sum()) * 1e6
merged_all_500[f'RNA_rep_comb_neuron_normalized500'] = (merged_all_500[f'RNA_rep_comb_neuron'] / merged_all_500[f'RNA_rep_comb_neuron'].sum()) * 1e6

merged_all_500["ratio_normalized500_chondro"] = np.log2(merged_all_500[f'RNA_rep_comb_chondro_normalized500']/merged_all_500[f'DNA_rep_comb_chondro_normalized500'])
merged_all_500["ratio_normalized500_neuron"] = np.log2(merged_all_500[f'RNA_rep_comb_neuron_normalized500']/merged_all_500[f'DNA_rep_comb_neuron_normalized500'])

merged_all_500_no_na = merged_all_500.dropna(subset=['ratio_normalized500_chondro', 'ratio_normalized500_neuron'], how='any')

plt.clf()
values = np.vstack([merged_all_500_no_na[f"ratio_normalized500_neuron"], merged_all_500_no_na[f"ratio_normalized500_chondro"]])
kernel = stats.gaussian_kde(values)(values)
sns.scatterplot(data = merged_all_500_no_na, x = f"ratio_normalized500_neuron", y = f"ratio_normalized500_chondro", s=1, c=kernel)
plt.axes().set_aspect('equal')
plt.savefig(f'./additional/chondro_neuron_comparison/ratio_chondro_vs_neurons_scatter_500dna_renormalized_{library}{adaptor}.pdf')
plt.savefig(f'./additional/chondro_neuron_comparison/ratio_chondro_vs_neurons_scatter_500dna_renormalized_{library}{adaptor}.png', dpi=330)