    # this script takes the comb df files, merges it with a dataframe that contains sequence and gc content and then runs analysis related to activity, gc content and coverage

import pandas as pd
from Bio import SeqIO
import sys
import os
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

activity_folder_path = f'./{cells}/{library}{adaptor}/output/activity_after_filter'
comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv' # N.M updated comb_df.csv to comb_df_adjusted_fdr.csv
fasta_file = f"./oligo_fasta/{library}{adaptor}.fasta"

comb_df = pd.read_csv(comb_df_path)

# Initialize empty lists to store modified identifiers and sequences
identifiers = []
sequences = []
gc_contents = []

# Parse the FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    identifier = record.id  # Remove the first character
    identifiers.append(identifier)
    sequence = str(record.seq)
    sequences.append(sequence)
    
    # Calculate GC content
    gc_count = sequence.count('g') + sequence.count('c')+sequence.count('G') + sequence.count('C')
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    gc_contents.append(gc_content)

# Create a DataFrame
gc_df = pd.DataFrame({'oligo': identifiers, 'Sequence': sequences, 'GC_Content': gc_contents})
print(gc_df.head().to_string())

# merge this to comb_df
merged_gc_comb = comb_df.merge(gc_df, on = "oligo", how = "left")

# turn dna count column into log scale to use for plotting (hue)
#merged_gc_comb_DNA5 = merged_gc_comb[merged_gc_comb["DNA_rep_comb"] >= 5 ]
merged_gc_comb_DNA5 = merged_gc_comb
merged_gc_comb_DNA5.loc[:,'DNA_rep_comb_log10'] = np.log10(merged_gc_comb_DNA5['DNA_rep_comb'])
merged_gc_comb_DNA5.loc[:,'DNA_rep_comb_clipped_500'] = merged_gc_comb_DNA5['DNA_rep_comb'].clip(upper=500, inplace=False)
merged_gc_comb_DNA5.loc[:, 'GC_Content_clipped_25_60'] = merged_gc_comb_DNA5['GC_Content'].clip(lower=25, upper=60, inplace=False)

# # plot correlation and color according to dna counts
# plt.clf()
# sns.scatterplot(data = merged_gc_comb_DNA5[(merged_gc_comb_DNA5[f"DNA_rep1"] >= 5)& (merged_gc_comb_DNA5[f"DNA_rep2"] >= 5)], x = f"ratio_log_rep1", y = f"ratio_log_rep2", hue="DNA_rep_comb_log10", s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - >=5bcs, color - DNA counts')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_DNA_counts_log10.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_DNA_counts_log10.png', dpi=330)

# # plot correlation and color according to dna counts
# plt.clf()
# sns.scatterplot(data = merged_gc_comb_DNA5[(merged_gc_comb_DNA5[f"DNA_rep1"] >= 5)& (merged_gc_comb_DNA5[f"DNA_rep2"] >= 5)], x = f"ratio_log_rep1", y = f"ratio_log_rep2", hue="DNA_rep_comb_clipped_500", s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - >=5bcs, color - DNA counts')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_DNA_counts_clipped500.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_DNA_counts_clipped500.png', dpi=330)

# # plot correlation and color according to gc_content
# plt.clf()
# sns.scatterplot(data = merged_gc_comb_DNA5[(merged_gc_comb_DNA5[f"DNA_rep1"] >= 5)& (merged_gc_comb_DNA5[f"DNA_rep2"] >= 5)], x = f"ratio_log_rep1", y = f"ratio_log_rep2", hue="GC_Content", s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - >=5bcs, color - GC')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_GC.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_GC.png', dpi=330)

# # plot correlation and color according to gc_content
# plt.clf()
# sns.scatterplot(data = merged_gc_comb_DNA5[(merged_gc_comb_DNA5[f"DNA_rep1"] >= 5)& (merged_gc_comb_DNA5[f"DNA_rep2"] >= 5)], x = f"ratio_log_rep1", y = f"ratio_log_rep2", hue="GC_Content_clipped_25_60", s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - >=5bcs, color - GC')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_GC_clipped2560.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_5_DNA_rna_dna_GC_clipped2560.png', dpi=330)


# take oligos with specifc DNA count and compare high vs low gc
bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]
labels = [
    '0-5', '5-10', '10-15', '15-20', '20-25', '25-30', '30-35', '35-40',
    '40-45', '45-50', '50-55', '55-60', '60-65', '65-70', '70-75', '75-80',
    '80-85', '85-90', '90-95', '95-100'
]


merged_gc_comb_DNA5.loc[:,'GC_Content_label'] = pd.cut(merged_gc_comb_DNA5['GC_Content'], bins=bins, labels=labels, right=True)
merged_gc_comb_DNA5.to_csv(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_df.csv', header=True, index = False, columns=["oligo", "alpha", "GC_Content", "GC_Content_label","GC_Content_clipped_25_60", "DNA_rep_comb", "DNA_rep_comb_clipped_500", "DNA_rep_comb_log10", "ratio_log_rep_comb"])

################################ RNA DNA log ratio plot##########################

f, (ax_hist,ax_bar) = plt.subplots(2, sharex=True,gridspec_kw={"height_ratios": (.5, .5)})

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["ratio_log_rep_comb"],
    showfliers=False,
    color="dodgerblue",
    ax=ax_hist)

sns.histplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    color="dodgerblue",
    ax=ax_bar)

ax_bar.tick_params(axis='x', labelsize=6,rotation=45)
ax_hist.set_ylabel("Reads per oligo")
ax_bar.set_ylabel("Number of oligos")
ax_hist.set_title(f"Reads per oligo by oligo gc content - {library}{adaptor}")
ax_bar.bar_label(ax_bar.containers[0],fontsize=5)

f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_ratio.pdf')
f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_ratio.png', dpi=330)

################################ combined plot ##########################

f, (ax_DNA, ax_RNA ,ax_hist) = plt.subplots(3, sharex=True,gridspec_kw={"height_ratios": (.33, .33, .34)})

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["DNA_rep_comb"],
    showfliers=False,
    color="dodgerblue",
    ax=ax_DNA)

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["RNA_rep_comb"],
    showfliers=False,
    color="dodgerblue",
    ax=ax_RNA)
       
    
sns.histplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    color="dodgerblue",
    ax=ax_hist)

ax_hist.tick_params(axis='x', labelsize=6,rotation=45)
ax_DNA.set_ylabel("DNA reads")
ax_RNA.set_ylabel("RNA reads")
ax_hist.set_ylabel("Number of oligos")
#ax_hist.set_title(f"Reads per oligo by oligo gc content - {library}{adaptor}")
ax_hist.bar_label(ax_bar.containers[0],fontsize=5)

f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_final.pdf')
f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_final.png', dpi=330)


################################ log combined plot ##########################

f, (ax_DNA, ax_RNA ,ax_hist) = plt.subplots(3, sharex=True,gridspec_kw={"height_ratios": (.33, .33, .34)})

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=np.log10(merged_gc_comb_DNA5["DNA_rep_comb"]),
    showfliers=False,
    color="dodgerblue",
    ax=ax_DNA)

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=np.log10(merged_gc_comb_DNA5["RNA_rep_comb"]),
    showfliers=False,
    color="dodgerblue",
    ax=ax_RNA)
       
    
sns.histplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    color="dodgerblue",
    ax=ax_hist)

ax_hist.tick_params(axis='x', labelsize=6,rotation=45)
ax_DNA.set_ylabel("log DNA reads")
ax_RNA.set_ylabel("log RNA reads")
ax_hist.set_ylabel("Number of oligos")
#ax_hist.set_title(f"Reads per oligo by oligo gc content - {library}{adaptor}")
ax_hist.bar_label(ax_bar.containers[0],fontsize=5)

f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_log_final.pdf')
f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_log_final.png', dpi=330)



################################  combined plot with normalized counts wo outliers ##########################


DNA_counts = merged_gc_comb_DNA5["DNA_rep_comb"].sum()
RNA_counts = merged_gc_comb_DNA5["RNA_rep_comb"].sum()

print(f"DNA counts: {DNA_counts}")
print(f"RNA counts: {RNA_counts}")

merged_gc_comb_DNA5["norm_DNA"] = ((merged_gc_comb_DNA5["DNA_rep_comb"]+1)*1000000)/DNA_counts
merged_gc_comb_DNA5["norm_RNA"] = ((merged_gc_comb_DNA5["RNA_rep_comb"]+1)*1000000)/RNA_counts

f, (ax_DNA, ax_RNA ,ax_hist) = plt.subplots(3, sharex=True,gridspec_kw={"height_ratios": (.33, .33, .34)})

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["norm_DNA"],
    showfliers=False,
    color="dodgerblue",
    ax=ax_DNA)

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["norm_RNA"],
    showfliers=False,
    color="dodgerblue",
    ax=ax_RNA)
       
    
sns.histplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    color="dodgerblue",
    ax=ax_hist)

ax_hist.tick_params(axis='x', labelsize=6,rotation=45)
ax_DNA.set_ylabel("DNA reads (CPM)")
ax_RNA.set_ylabel("RNA reads (CPM)")
ax_hist.set_ylabel("Number of oligos")
#ax_hist.set_title(f"Reads per oligo by oligo gc content - {library}{adaptor}")
ax_hist.bar_label(ax_bar.containers[0],fontsize=5)

f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_normalized_wo_outliers_final.pdf')
f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_normalized_wo_outliers_final.png', dpi=330)

################################  combined plot with normalized counts##########################


DNA_counts = merged_gc_comb_DNA5["DNA_rep_comb"].sum()
RNA_counts = merged_gc_comb_DNA5["RNA_rep_comb"].sum()

print(f"DNA counts: {DNA_counts}")
print(f"RNA counts: {RNA_counts}")

merged_gc_comb_DNA5["norm_DNA"] = ((merged_gc_comb_DNA5["DNA_rep_comb"]+1)*1000000)/DNA_counts
merged_gc_comb_DNA5["norm_RNA"] = ((merged_gc_comb_DNA5["RNA_rep_comb"]+1)*1000000)/RNA_counts

f, (ax_DNA, ax_RNA ,ax_hist) = plt.subplots(3, sharex=True,gridspec_kw={"height_ratios": (.33, .33, .34)})

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["norm_DNA"],
    showfliers=True,
    color="dodgerblue",
    ax=ax_DNA)

sns.boxplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    y=merged_gc_comb_DNA5["norm_RNA"],
    showfliers=True,
    color="dodgerblue",
    ax=ax_RNA)
       
    
sns.histplot(
    x=merged_gc_comb_DNA5["GC_Content_label"],
    color="dodgerblue",
    ax=ax_hist)

ax_hist.tick_params(axis='x', labelsize=6,rotation=45)
ax_DNA.set_ylabel("DNA reads (CPM)")
ax_RNA.set_ylabel("RNA reads (CPM)")
ax_hist.set_ylabel("Number of oligos")
#ax_hist.set_title(f"Reads per oligo by oligo gc content - {library}{adaptor}")
ax_hist.bar_label(ax_bar.containers[0],fontsize=5)

f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_normalized_final.pdf')
f.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_normalized_final.png', dpi=330)



################################ DNA counts boxplot plot ##########################

#merged_gc_comb_DNA5_9 = merged_gc_comb_DNA5[merged_gc_comb_DNA5["DNA_rep_comb"] <10].copy()
plt.clf()
ax = sns.boxplot(x=merged_gc_comb_DNA5["GC_Content_label"], y=merged_gc_comb_DNA5["DNA_rep_comb"],showfliers=False)

# Calculate number of obs per group & median to position labels
medians = merged_gc_comb_DNA5.groupby(['GC_Content_label'])['DNA_rep_comb'].median().values
nobs = merged_gc_comb_DNA5['GC_Content_label'].value_counts(sort=False).values
print(nobs)
nobs = [str(x) for x in nobs.tolist()]
print(nobs)
nobs = ["n: " + i for i in nobs]
print(nobs)

# Add text to the figure
pos = range(len(nobs))
print(pos)
for tick, label in zip(pos, ax.get_xticklabels()):
    print(tick)
    print(label)
    ax.text(pos[tick], medians[tick] + 0.03, nobs[tick],
            horizontalalignment='center',
            size='xx-small',
            color='y',
            weight='semibold')

plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_boxplot_gc_DNA_counts.pdf')
plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_boxplot_gc_DNA_counts.png', dpi=330)

################################ Ratio boxplot ##########################

merged_gc_comb_DNA5_9 = merged_gc_comb_DNA5[merged_gc_comb_DNA5["DNA_rep_comb"] <10].copy()
plt.clf()
ax = sns.boxplot(x=merged_gc_comb_DNA5["GC_Content_label"], y=merged_gc_comb_DNA5["ratio_log_rep_comb"],showfliers=False)

# Calculate number of obs per group & median to position labels
medians = merged_gc_comb_DNA5.groupby(['GC_Content_label'])['ratio_log_rep_comb'].median().values
nobs = merged_gc_comb_DNA5['GC_Content_label'].value_counts(sort=False).values
print(nobs)
nobs = [str(x) for x in nobs.tolist()]
print(nobs)
nobs = ["n: " + i for i in nobs]
print(nobs)

# Add text to the figure
pos = range(len(nobs))
print(pos)
for tick, label in zip(pos, ax.get_xticklabels()):
    print(tick)
    print(label)
    ax.text(pos[tick], medians[tick] + 0.03, nobs[tick],
            horizontalalignment='center',
            size='xx-small',
            color='y',
            weight='semibold')

plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_boxplot_gc_activity_ratio.pdf')
plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_boxplot_gc_activity_ratio.png', dpi=330)

################################ RNA counts boxplot ##########################

merged_gc_comb_DNA5_9 = merged_gc_comb_DNA5[merged_gc_comb_DNA5["DNA_rep_comb"] <10].copy()
plt.clf()
ax = sns.boxplot(x=merged_gc_comb_DNA5["GC_Content_label"], y=merged_gc_comb_DNA5["RNA_rep_comb"],showfliers=False)

# Calculate number of obs per group & median to position labels
medians = merged_gc_comb_DNA5.groupby(['GC_Content_label'])['RNA_rep_comb'].median().values
nobs = merged_gc_comb_DNA5['GC_Content_label'].value_counts(sort=False).values
print(nobs)
nobs = [str(x) for x in nobs.tolist()]
print(nobs)
nobs = ["n: " + i for i in nobs]
print(nobs)

# Add text to the figure
pos = range(len(nobs))
print(pos)
for tick, label in zip(pos, ax.get_xticklabels()):
    print(tick)
    print(label)
    ax.text(pos[tick], medians[tick] + 0.03, nobs[tick],
            horizontalalignment='center',
            size='xx-small',
            color='y',
            weight='semibold')

plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_boxplot_gc_activity_RNA.pdf')
plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}/{library}{adaptor}_boxplot_gc_activity_RNA.png', dpi=330)


# merged_gc_comb_DNA50_70 = merged_gc_comb_DNA5[(merged_gc_comb_DNA5["DNA_rep_comb"] <=70) & (merged_gc_comb_DNA5["DNA_rep_comb"] >=50)].copy()
# plt.clf()
# ax = sns.violinplot(x=merged_gc_comb_DNA50_70["GC_Content_label"], y=merged_gc_comb_DNA50_70["ratio_log_rep_comb"])

# # Calculate number of obs per group & median to position labels
# medians = merged_gc_comb_DNA50_70.groupby(['GC_Content_label'])['ratio_log_rep_comb'].median().values
# nobs = merged_gc_comb_DNA50_70['GC_Content_label'].value_counts(sort=False).values
# nobs = [str(x) for x in nobs.tolist()]
# nobs = ["n: " + i for i in nobs]
 
# # Add text to the figure
# pos = range(len(nobs))
# for tick, label in zip(pos, ax.get_xticklabels()):
    # ax.text(pos[tick], medians[tick] + 0.03, nobs[tick],
            # horizontalalignment='center',
            # size='xx-small',
            # color='w',
            # weight='semibold')

# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_violin_gc_activity_ratio_dna50_70.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_violin_gc_activity_ratio_dna50_70.png', dpi=330)


# create a measure of noise (distance from diagonal x=y) and plot noise in relation to other features (eg gc, dna, alpha ... )
merged_gc_comb_DNA5["noise"]= abs(merged_gc_comb_DNA5["ratio_log_rep1"] - merged_gc_comb_DNA5["ratio_log_rep2"])/np.sqrt(2)
# merged_gc_comb_DNA5.to_csv(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_noise_df.csv', header=True, index = False, columns=["oligo", "alpha", "GC_Content", "DNA_rep_comb", "DNA_rep1", "DNA_rep2", "DNA_rep3", "noise", "ratio_log_rep_comb"])

# plt.clf()
# plt.hist(merged_gc_comb_DNA5["noise"], bins=100)
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_noise_hist_100.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_noise_hist_100.png', dpi=330)

merged_gc_comb_DNA5_no_na = merged_gc_comb_DNA5.dropna(subset=['noise'], inplace=False)

# plt.clf()
# values = np.vstack([merged_gc_comb_DNA5_no_na["GC_Content"], merged_gc_comb_DNA5_no_na["noise"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_gc_comb_DNA5_no_na, x = f"GC_Content", y = f"noise", c=kernel, s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation noise GC')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_GC_color.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_GC_color.png', dpi=330)

# corr = merged_gc_comb_DNA5_no_na["GC_Content"].corr(merged_gc_comb_DNA5_no_na["noise"])
# print("gc-noise", corr)

# plt.clf()
# values = np.vstack([merged_gc_comb_DNA5_no_na[f"DNA_rep_comb"], merged_gc_comb_DNA5_no_na[f"noise"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_gc_comb_DNA5_no_na, x = f"DNA_rep_comb", y = f"noise", c=kernel, s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation noise DNA')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_DNA_color.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_DNA_color.png', dpi=330)

# corr = merged_gc_comb_DNA5_no_na["DNA_rep_comb"].corr(merged_gc_comb_DNA5_no_na["noise"])
# print("DNA-noise", corr)

# plt.clf()
# values = np.vstack([merged_gc_comb_DNA5_no_na[f"ratio_log_rep_comb"], merged_gc_comb_DNA5_no_na[f"noise"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_gc_comb_DNA5_no_na, x = f"ratio_log_rep_comb", y = f"noise", c=kernel, s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation noise activity')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_activity_color.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_activity_color.png', dpi=330)

# plt.clf()
# sns.scatterplot(data = merged_gc_comb_DNA5_no_na, x = f"ratio_log_rep_comb", y = f"noise", hue="DNA_rep_comb_log10", s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation noise activity')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_activity_color_coverage.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_noise_activity_color_coverage.png', dpi=330)


# corr = merged_gc_comb_DNA5_no_na["ratio_log_rep_comb"].corr(merged_gc_comb_DNA5_no_na["noise"])
# print("activity-noise", corr)

# plt.clf()
# values = np.vstack([merged_gc_comb_DNA5_no_na[f"GC_Content"], merged_gc_comb_DNA5_no_na[f"DNA_rep_comb"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = merged_gc_comb_DNA5_no_na, x = f"GC_Content", y = f"DNA_rep_comb", c=kernel, s=1)
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation DNA GC')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_GC_DNA_color.pdf')
# plt.savefig(f'./additional/GC_activity_coverage/{cells}/{library}{adaptor}_scatter_GC_DNA_color.png', dpi=330)

# corr = merged_gc_comb_DNA5_no_na["GC_Content"].corr(merged_gc_comb_DNA5_no_na["DNA_rep_comb"])
# print("gc-dna", corr)