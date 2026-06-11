import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from scipy import stats
import matplotlib.colors as mcolors


library = sys.argv[1]
adaptor = sys.argv[2]

comb_df_path_neurons = f'./neurons/{library}{adaptor}/output/activity_after_filter/comb_df_adjusted_fdr.csv'
comb_df_path_chondro = f'./chondrocytes/{library}{adaptor}/output/activity_after_filter/comb_df_adjusted_fdr.csv'

neuron_df = pd.read_csv(comb_df_path_neurons)
chondro_df = pd.read_csv(comb_df_path_chondro)
umi_chondro = pd.read_csv(f'./chondrocytes/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t', header=0)
umi_neuron = pd.read_csv(f'./neurons/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t', header=0)

merged_all = neuron_df.merge(chondro_df, on="oligo", suffixes=('_neuron', '_chondro'))

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



for rep in ["RNA_rep1","RNA_rep2","RNA_rep3","RNA_rep_comb"]:
    neuron_df[f'{rep}_normalized'] = (neuron_df[f'{rep}'] / neuron_df[f'{rep}'].sum()) * 1e6
    chondro_df[f'{rep}_normalized'] = (chondro_df[f'{rep}'] / chondro_df[f'{rep}'].sum()) * 1e6
for rep in ["DNA_rep1","DNA_rep2","DNA_rep3","DNA_rep_comb"]:
    neuron_df[f'{rep}_normalized'] = (neuron_df[f'{rep}'] / neuron_df[f'{rep}'].sum()) * 1e6
    chondro_df[f'{rep}_normalized'] = (chondro_df[f'{rep}'] / chondro_df[f'{rep}'].sum()) * 1e6

merged_all_with_normalized = neuron_df.merge(chondro_df, on="oligo", suffixes=('_neuron', '_chondro'))

# histogram of rna counts per oligo comparing neurons and chondrocytes + scatter (counts normalized to library size)
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
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_RNA_distribution_cpm_clipped125_{library}{adaptor}.png', transparent=True, dpi=330)
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_RNA_distribution_cpm_clipped125_{library}{adaptor}.pdf')

# merged_all_with_normalized_no_na_rna = merged_all_with_normalized.dropna(subset=['RNA_rep_comb_normalized_neuron', 'RNA_rep_comb_normalized_chondro'], how='any')

# plt.clf()
# values = np.vstack([merged_all_with_normalized_no_na_rna[f"RNA_rep_comb_normalized_neuron"], merged_all_with_normalized_no_na_rna[f"RNA_rep_comb_normalized_chondro"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_all_with_normalized_no_na_rna, x = f"RNA_rep_comb_normalized_chondro", y = f"RNA_rep_comb_normalized_neuron", s=1, c=kernel)
# plt.axes().set_aspect('equal')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_rna_chondro_vs_neurons_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_rna_chondro_vs_neurons_scatter_{library}{adaptor}.png', dpi=330)

# # histogram of dna counts per oligo comparing neurons and chondrocytes + scatter (counts normalized to library size)


# plt.clf()
# fig, axes = plt.subplots(2,3, sharey=True, sharex=True)
# for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    # curr_df_neurons = neuron_df.copy()
    # curr_df_neurons[f"DNA_{rep}_normalized"] = curr_df_neurons[f"DNA_{rep}_normalized"].clip(upper=125)
    # curr_df_chondro = chondro_df.copy()
    # curr_df_chondro[f"DNA_{rep}_normalized"] = curr_df_chondro[f"DNA_{rep}_normalized"].clip(upper=125)
    # sns.histplot(ax=axes[0, n], data=curr_df_neurons, x=f"DNA_{rep}_normalized", edgecolor='none')
    # sns.histplot(ax=axes[1, n], data=curr_df_chondro, x=f"DNA_{rep}_normalized", edgecolor='none')
# plt.suptitle("Distribution DNA (cpm) - clipped to 125")
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_DNA_distribution_cpm_clipped125_{library}{adaptor}.png', transparent=True, dpi=330)
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_DNA_distribution_cpm_clipped125_{library}{adaptor}.pdf')

# merged_all_with_normalized_no_na_DNA = merged_all_with_normalized.dropna(subset=['DNA_rep_comb_normalized_neuron', 'DNA_rep_comb_normalized_chondro'], how='any')

# plt.clf()
# values = np.vstack([merged_all_with_normalized_no_na_DNA[f"DNA_rep_comb_normalized_neuron"], merged_all_with_normalized_no_na_DNA[f"DNA_rep_comb_normalized_chondro"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_all_with_normalized_no_na_DNA, x = f"DNA_rep_comb_normalized_chondro", y = f"DNA_rep_comb_normalized_neuron", s=1, c=kernel)
# plt.axes().set_aspect('equal')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_DNA_chondro_vs_neurons_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/oligo_DNA_chondro_vs_neurons_scatter_{library}{adaptor}.png', dpi=330)


# # correlations between dna and rna counts per oligo (for chondrocytes and neurons)
# neuron_df_no_na = neuron_df.dropna(subset=['DNA_rep1_normalized', 'RNA_rep1_normalized'], how='any')
# chondro_df_no_na = chondro_df.dropna(subset=['DNA_rep1_normalized', 'RNA_rep1_normalized'], how='any')
# max_chondro_rna = max(chondro_df_no_na["RNA_rep1_normalized"]) +50
# print(max_chondro_rna)
# plt.clf()
# fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)
# values_neurons = np.vstack([neuron_df_no_na[f"DNA_rep1_normalized"], neuron_df_no_na[f"RNA_rep1_normalized"]])
# kernel_neurons = stats.gaussian_kde(values_neurons)(values_neurons)
# values_chondro = np.vstack([chondro_df_no_na[f"DNA_rep1_normalized"], chondro_df_no_na[f"RNA_rep1_normalized"]])
# kernel_chondro = stats.gaussian_kde(values_chondro)(values_chondro)
# axs[0].set_xlim(0.0001, max_chondro_rna)
# axs[1].set_xlim(0.0001, max_chondro_rna)
# axs[0].set_ylim(0.0001, max_chondro_rna)
# axs[1].set_ylim(0.0001, max_chondro_rna)
# axs[0].set_yscale('log', base=2)
# axs[1].set_yscale('log', base=2)
# axs[0].set_xscale('log', base=2)
# axs[1].set_xscale('log', base=2)
# sns.scatterplot(ax=axs[0], data=neuron_df_no_na, x=f"DNA_rep1_normalized", y="RNA_rep1_normalized", s=1, linewidth=0, c=kernel_neurons)
# sns.scatterplot(ax=axs[1], data=chondro_df_no_na, x=f"DNA_rep1_normalized", y="RNA_rep1_normalized", s=1, linewidth=0, c=kernel_chondro)
# axs[0].set_title('Neurons')
# axs[1].set_title('Chondrocytes')
# axs[0].set_aspect('equal')
# axs[1].set_aspect('equal')
# plt.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/rna_vs_dna_oligo_counts_scatter_{library}{adaptor}_color_log2.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/rna_vs_dna_oligo_counts_scatter_{library}{adaptor}_color_log2.png', dpi=330)

# corr_neurons = neuron_df_no_na['DNA_rep1_normalized'].corr(neuron_df_no_na["RNA_rep1_normalized"])
# corr_chondro = chondro_df_no_na['DNA_rep1_normalized'].corr(chondro_df_no_na["RNA_rep1_normalized"])

# print(f'oligo - corr dna(chondrocytes)- rna(chondrocytes) {corr_chondro}')
# print(f'oligo - corr dna(neurons)- rna(neurons) {corr_neurons}')

# # correlations between dna and rna counts per barcode (for chondrocytes and neurons)
# max_chondro_rna_bc = max(umi_chondro["rna_rep1_normalized"])+1
# print(max_chondro_rna_bc)
# plt.clf()
# fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)
# sns.scatterplot(ax=axs[0], data=umi_neuron, x=f"dna_rep1_normalized", y="rna_rep1_normalized", s=1, linewidth=0)
# sns.scatterplot(ax=axs[1], data=umi_chondro, x=f"dna_rep1_normalized", y="rna_rep1_normalized", s=1, linewidth=0)
# axs[0].set_title('Neurons')
# axs[1].set_title('Chondrocytes')
# axs[0].set_aspect('equal')
# axs[1].set_aspect('equal')
# axs[0].set_xlim(right=max_chondro_rna_bc)
# axs[1].set_xlim(right=max_chondro_rna_bc)
# axs[0].set_ylim(top=max_chondro_rna_bc)
# axs[1].set_ylim(top=max_chondro_rna_bc)
# plt.tight_layout()
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/rna_vs_dna_barcode_counts_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/rna_vs_dna_barcode_counts_scatter_{library}{adaptor}.png', dpi=330)

# corr_neurons = umi_neuron['dna_rep1_normalized'].corr(umi_neuron["rna_rep1_normalized"])
# corr_chondro = umi_chondro['dna_rep1_normalized'].corr(umi_chondro["rna_rep1_normalized"])

# print(f'bc - corr dna(chondrocytes)- rna(chondrocytes) {corr_chondro}')
# print(f'bc - corr dna(neurons)- rna(neurons) {corr_neurons}')

# # per barcode: dna-dna, rna-rna (chondro vs neurons)
# umi_merged = umi_neuron.merge(umi_chondro, on="oligo_bc", suffixes=('_neuron', '_chondro'))
# plt.clf()
# sns.scatterplot(data=umi_merged, x=f"dna_rep1_normalized_chondro", y="dna_rep1_normalized_neuron", s=1, linewidth=0)
# plt.axes().set_aspect('equal')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_dna_counts_neurons_vs_chondro_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_dna_counts_neurons_vs_chondro_scatter_{library}{adaptor}.png', dpi=330)

# plt.clf()
# sns.scatterplot(data=umi_merged, x=f"rna_rep1_normalized_chondro", y="rna_rep1_normalized_neuron", s=1, linewidth=0)
# plt.axes().set_aspect('equal')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_rna_counts_neurons_vs_chondro_scatter_{library}{adaptor}.pdf')
# plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_rna_counts_neurons_vs_chondro_scatter_{library}{adaptor}.png', dpi=330)

# d_c_d_n = umi_merged['dna_rep1_normalized_chondro'].corr(umi_merged["dna_rep1_normalized_neuron"])
# r_c_r_n = umi_merged['rna_rep1_normalized_chondro'].corr(umi_merged["rna_rep1_normalized_neuron"])

# print(f'corr dna(chondrocytes)- dna(neurons) {d_c_d_n}')
# print(f'corr rna(chondrocytes)- rna(neurons) {r_c_r_n}')

# heat map correlations between dna and rna counts per barcode (for chondrocytes and neurons)
max_chondro_rna_bc = max(umi_chondro["rna_rep1_normalized"])+1
print(max_chondro_rna_bc)
plt.clf()
fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)
hist1 = axs[0].hist2d(umi_neuron["dna_rep1_normalized"], umi_neuron["rna_rep1_normalized"], bins=100, cmap = 'Reds', range=[[0,max_chondro_rna_bc],[0,max_chondro_rna_bc]], norm=mcolors.LogNorm())
hist2 = axs[1].hist2d(umi_chondro["dna_rep1_normalized"], umi_chondro["rna_rep1_normalized"], bins=100, cmap = 'Reds', range=[[0,max_chondro_rna_bc],[0,max_chondro_rna_bc]], norm=mcolors.LogNorm())
axs[0].set_title('Neurons')
axs[1].set_title('Chondrocytes')
axs[0].set_aspect('equal')
axs[1].set_aspect('equal')
axs[0].set_xlim(right=max_chondro_rna_bc)
axs[1].set_xlim(right=max_chondro_rna_bc)
axs[0].set_ylim(top=max_chondro_rna_bc)
axs[1].set_ylim(top=max_chondro_rna_bc)
axs[0].set_xlabel("dna cpm rep1")
axs[1].set_xlabel("dna cpm rep1")
axs[0].set_ylabel("rna cpm rep1")
axs[1].set_ylabel("rna cpm rep1")
cbar1 = fig.colorbar(hist1[3], ax=axs[0])
cbar2 = fig.colorbar(hist2[3], ax=axs[1])
plt.tight_layout()
plt.savefig(f'./additional/chondro_neuron_comparison_thesis/rna_vs_dna_barcode_counts_heatmap_{library}{adaptor}.pdf')
plt.savefig(f'./additional/chondro_neuron_comparison_thesis/rna_vs_dna_barcode_counts_heatmap_{library}{adaptor}.png', dpi=330)

# # heatmap per barcode: dna-dna, rna-rna (chondro vs neurons)
umi_merged = umi_neuron.merge(umi_chondro, on="oligo_bc", suffixes=('_neuron', '_chondro'))
plt.clf()
plt.hist2d(umi_merged["dna_rep1_normalized_chondro"], umi_merged["dna_rep1_normalized_neuron"], bins=100, cmap = 'Reds', norm=mcolors.LogNorm())
plt.axes().set_aspect('equal')
plt.xlabel('chondro - dna counts per barcode, rep1')
plt.ylabel('neuron - dna counts per barcode, rep1')
plt.tight_layout()
cb = plt.colorbar()
cb.set_label('counts in bin')
plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_dna_counts_neurons_vs_chondro_heatmap_{library}{adaptor}.pdf')
plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_dna_counts_neurons_vs_chondro_heatmap_{library}{adaptor}.png', dpi=330)

plt.clf()
plt.hist2d(umi_merged["rna_rep1_normalized_chondro"], umi_merged["rna_rep1_normalized_neuron"], bins=100, cmap = 'Reds', norm=mcolors.LogNorm())
plt.axes().set_aspect('equal')
plt.xlabel('chondro - Rna counts per barcode, rep1')
plt.ylabel('neuron - Rna counts per barcode, rep1')
plt.tight_layout()
cb = plt.colorbar()
cb.set_label('counts in bin')
plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_rna_counts_neurons_vs_chondro_heatmap_{library}{adaptor}.pdf')
plt.savefig(f'./additional/chondro_neuron_comparison_thesis/barcode_rna_counts_neurons_vs_chondro_heatmap_{library}{adaptor}.png', dpi=330)