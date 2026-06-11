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
cells = sys.argv[3]
percent = sys.argv[4]


#Obtaining the read counts with soft coding, Nadav 06-10-2024
UMI_all_reps = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t')
UMI_all_reps.set_index('oligo_bc', inplace=True)
rep1DNA=UMI_all_reps[f"DNA_{cells}_rep1"].sum()
rep1RNA=UMI_all_reps[f"RNA_{cells}_rep1"].sum()
rep2DNA=UMI_all_reps[f"DNA_{cells}_rep2"].sum()
rep2RNA=UMI_all_reps[f"RNA_{cells}_rep2"].sum()
rep3DNA=UMI_all_reps[f"DNA_{cells}_rep3"].sum()
rep3RNA=UMI_all_reps[f"RNA_{cells}_rep3"].sum()

print("rep1DNA:" ,rep1DNA)
print("rep1RNA:" ,rep1RNA)
print("rep2DNA:" ,rep2DNA)
print("rep2RNA:" ,rep2RNA)
print("rep3DNA:" ,rep3DNA)
print("rep3RNA:" ,rep3RNA)

DNA_reads = [rep1DNA, rep2DNA, rep3DNA]
RNA_reads = [rep1RNA, rep2RNA, rep3RNA]


# access replicate files
rep1 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep1.txt', sep='\t', header=0, names=["oligo_bc", "RNA", "DNA"])
rep2 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep2.txt', sep='\t', header=0, names=["oligo_bc", "RNA", "DNA"])
rep3 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep3.txt', sep='\t', header=0, names=["oligo_bc", "RNA", "DNA"])

# take out barcodes that have RNA and DNA 0
rep1 = rep1[(rep1['DNA'] > 0) | (rep1['RNA'] > 0)]
rep2 = rep2[(rep2['DNA'] > 0) | (rep2['RNA'] > 0)]
rep3 = rep3[(rep3['DNA'] > 0) | (rep3['RNA'] > 0)]


# group oligos and sum RNA, DNA counts. Also counts how many barcodes we have per oligo. Caluclate ratio and perform log2 on ratio
rep_list = []

for i, rep in enumerate([rep1, rep2, rep3]):
    rep[['oligo','bc']] = rep.oligo_bc.str.rsplit("_", n=1, expand=True)
    grouped_df = rep.groupby('oligo')
    final = grouped_df[['RNA', 'DNA']].agg('sum')
    final["count"] = grouped_df.size()
    final["RNA_pseudo"] = final["RNA"] + 1
    final["DNA_pseudo"] = final["DNA"] + 1
    final["RNA_cpm"] = (final["RNA_pseudo"]*1000000)/RNA_reads[i]
    final["DNA_cpm"] = (final["DNA_pseudo"]*1000000)/DNA_reads[i]
    final["ratio"] = final["RNA_cpm"]/final["DNA_cpm"]
    final["ratio_log"] = np.log2(final["ratio"])
    final['oligo_column'] = final.index
    rep_list.append(final)


rep1_counts = rep_list[0]
rep2_counts = rep_list[1]
rep3_counts = rep_list[2]

# merge different replicates
rep1_2 = rep1_counts.merge(rep2_counts, how="inner", on="oligo_column")
rep2_3 = rep2_counts.merge(rep3_counts, how="inner", on="oligo_column")
rep3_1 = rep3_counts.merge(rep1_counts, how="inner", on="oligo_column")

print(type(rep1_2))
print(rep1_2.head().to_string())

# add labels for controls
rep1_2.loc[rep1_2['oligo_column'].str.contains('NegCtrl_not_active'), 'Label'] = 'Negative control'
rep1_2.loc[rep1_2['oligo_column'].str.contains(f'PosCtrl_osteoblast_active'), 'Label'] = 'Positive control'
rep2_3.loc[rep2_3['oligo_column'].str.contains('NegCtrl_not_active'), 'Label'] = 'Negative control'
rep2_3.loc[rep2_3['oligo_column'].str.contains(f'PosCtrl_osteoblast_active'), 'Label'] = 'Positive control'
rep3_1.loc[rep3_1['oligo_column'].str.contains('NegCtrl_not_active'), 'Label'] = 'Negative control'
rep3_1.loc[rep3_1['oligo_column'].str.contains(f'PosCtrl_osteoblast_active'), 'Label'] = 'Positive control'

rep1_2['Label'] = rep1_2['Label'].fillna(value='No control')
rep2_3['Label'] = rep2_3['Label'].fillna(value='No control')
rep3_1['Label'] = rep3_1['Label'].fillna(value='No control')


rep1_2['Label'] = pd.Categorical(rep1_2.Label, categories=['No control', 'Positive control', 'Negative control'], ordered=True)
rep1_2.sort_values('Label', inplace = True)
rep2_3['Label'] = pd.Categorical(rep2_3.Label, categories=['No control', 'Positive control', 'Negative control'], ordered=True)
rep2_3.sort_values('Label', inplace = True)
rep3_1['Label'] = pd.Categorical(rep3_1.Label, categories=['No control', 'Positive control', 'Negative control'], ordered=True)
rep3_1.sort_values('Label', inplace = True)

# filter to take only oligos with specific amount of barcodes
rep10bcs = []
rep5bcs = []
rep_all_bcs = []

for rep in [rep1_2, rep2_3, rep3_1]:
    # rep.replace([np.inf, -np.inf], np.nan, inplace=True) #can probably delete this line
    # rep = rep.dropna()
    rep10bcs.append(rep[(rep["count_x"] >= 10) & (rep["count_y"] >= 10)])
    rep5bcs.append(rep[(rep["count_x"] >= 5) & (rep["count_y"] >= 5)])
    rep_all_bcs.append(rep)


print(f"Length of rep10bcs for rep1_2: {len(rep10bcs[0])}")
print(f"Length of rep5bcs for rep1_2: {len(rep10bcs[1])}")
print(f"Length of all for rep1_2: {len(rep10bcs[2])}")

# calculate overall correlation between replicates >=10 bcs or >=5 bcs

corr_1_2_log = rep_all_bcs[0]["ratio_log_x"].corr(rep_all_bcs[0]["ratio_log_y"])
corr_2_3_log = rep_all_bcs[1]["ratio_log_x"].corr(rep_all_bcs[1]["ratio_log_y"])
corr_3_1_log = rep_all_bcs[2]["ratio_log_x"].corr(rep_all_bcs[2]["ratio_log_y"])

corr_1_2_10bcs_log = rep10bcs[0]["ratio_log_x"].corr(rep10bcs[0]["ratio_log_y"])
corr_2_3_10bcs_log = rep10bcs[1]["ratio_log_x"].corr(rep10bcs[1]["ratio_log_y"])
corr_3_1_10bcs_log = rep10bcs[2]["ratio_log_x"].corr(rep10bcs[2]["ratio_log_y"])

corr_1_2_5bcs_log = rep5bcs[0]["ratio_log_x"].corr(rep5bcs[0]["ratio_log_y"])
corr_2_3_5bcs_log = rep5bcs[1]["ratio_log_x"].corr(rep5bcs[1]["ratio_log_y"])
corr_3_1_5bcs_log = rep5bcs[2]["ratio_log_x"].corr(rep5bcs[2]["ratio_log_y"])

print(corr_1_2_log)
print(corr_2_3_log)
print(corr_3_1_log)

print(corr_1_2_10bcs_log)
print(corr_2_3_10bcs_log)
print(corr_3_1_10bcs_log)

print(corr_1_2_5bcs_log)
print(corr_2_3_5bcs_log)
print(corr_3_1_5bcs_log)


# calculate overall correlation between replicates >=10 bcs or >=5 bcs for oligos with log RNA/DNA above one # NM 17-11-2024 Added the following section
print("\nWith activity threshold (log):\n")

rep10bcs = []
rep5bcs = []
rep_all_bcs = []

for rep in [rep1_2, rep2_3, rep3_1]:

    filtered_rep10 = rep[(rep["count_x"] >= 10) & (rep["count_y"] >= 10) & (rep["ratio_log_x"] >= 1) & (rep["ratio_log_y"] >= 1)]
    rep10bcs.append(filtered_rep10)

    filtered_rep5 = rep[(rep["count_x"] >= 5) & (rep["count_y"] >= 5) & (rep["ratio_log_x"] >= 1) & (rep["ratio_log_y"] >= 1)]
    rep5bcs.append(filtered_rep5)

    filtered_rep_all = rep[rep["ratio_log_x"] >= 1 & (rep["ratio_log_y"] >= 1)]
    rep_all_bcs.append(filtered_rep_all)
    
    
print(f"Length of rep10bcs for rep1_2: {len(rep10bcs[0])}")
print(f"Length of rep5bcs for rep1_2: {len(rep10bcs[1])}")
print(f"Length of all for rep1_2: {len(rep10bcs[2])}")

corr_1_2_log = rep_all_bcs[0]["ratio_log_x"].corr(rep_all_bcs[0]["ratio_log_y"])
corr_2_3_log = rep_all_bcs[1]["ratio_log_x"].corr(rep_all_bcs[1]["ratio_log_y"])
corr_3_1_log = rep_all_bcs[2]["ratio_log_x"].corr(rep_all_bcs[2]["ratio_log_y"])

corr_1_2_10bcs_log = rep10bcs[0]["ratio_log_x"].corr(rep10bcs[0]["ratio_log_y"])
corr_2_3_10bcs_log = rep10bcs[1]["ratio_log_x"].corr(rep10bcs[1]["ratio_log_y"])
corr_3_1_10bcs_log = rep10bcs[2]["ratio_log_x"].corr(rep10bcs[2]["ratio_log_y"])

corr_1_2_5bcs_log = rep5bcs[0]["ratio_log_x"].corr(rep5bcs[0]["ratio_log_y"])
corr_2_3_5bcs_log = rep5bcs[1]["ratio_log_x"].corr(rep5bcs[1]["ratio_log_y"])
corr_3_1_5bcs_log = rep5bcs[2]["ratio_log_x"].corr(rep5bcs[2]["ratio_log_y"])

print(corr_1_2_log)
print(corr_2_3_log)
print(corr_3_1_log)

print(corr_1_2_10bcs_log)
print(corr_2_3_10bcs_log)
print(corr_3_1_10bcs_log)

print(corr_1_2_5bcs_log)
print(corr_2_3_5bcs_log)
print(corr_3_1_5bcs_log)



# calculate overall correlation between replicates >=10 bcs or >=5 bcs for oligos with log RNA/DNA above one # NM 17-11-2024 Added the following section
print("\nWith activity threshold:\n")

rep10bcs = []
rep5bcs = []
rep_all_bcs = []

for rep in [rep1_2, rep2_3, rep3_1]:

    filtered_rep10 = rep[(rep["count_x"] >= 10) & (rep["count_y"] >= 10) & (rep["ratio_log_x"] >= 1) & (rep["ratio_log_y"] >= 1)]
    rep10bcs.append(filtered_rep10)

    filtered_rep5 = rep[(rep["count_x"] >= 5) & (rep["count_y"] >= 5) & (rep["ratio_log_x"] >= 1) & (rep["ratio_log_y"] >= 1)]
    rep5bcs.append(filtered_rep5)

    filtered_rep_all = rep[rep["ratio_log_x"] >= 1 & (rep["ratio_log_y"] >= 1)]
    rep_all_bcs.append(filtered_rep_all)
    
    
print(f"Length of rep10bcs for rep1_2: {len(rep10bcs[0])}")
print(f"Length of rep5bcs for rep1_2: {len(rep10bcs[1])}")
print(f"Length of all for rep1_2: {len(rep10bcs[2])}")

corr_1_2 = rep_all_bcs[0]["ratio_x"].corr(rep_all_bcs[0]["ratio_y"])
corr_2_3 = rep_all_bcs[1]["ratio_x"].corr(rep_all_bcs[1]["ratio_y"])
corr_3_1 = rep_all_bcs[2]["ratio_x"].corr(rep_all_bcs[2]["ratio_y"])

corr_1_2_10bcs = rep10bcs[0]["ratio_x"].corr(rep10bcs[0]["ratio_y"])
corr_2_3_10bcs = rep10bcs[1]["ratio_x"].corr(rep10bcs[1]["ratio_y"])
corr_3_1_10bcs = rep10bcs[2]["ratio_x"].corr(rep10bcs[2]["ratio_y"])

corr_1_2_5bcs = rep5bcs[0]["ratio_x"].corr(rep5bcs[0]["ratio_y"])
corr_2_3_5bcs = rep5bcs[1]["ratio_x"].corr(rep5bcs[1]["ratio_y"])
corr_3_1_5bcs = rep5bcs[2]["ratio_x"].corr(rep5bcs[2]["ratio_y"])

print(corr_1_2)
print(corr_2_3)
print(corr_3_1)

print(corr_1_2_10bcs)
print(corr_2_3_10bcs)
print(corr_3_1_10bcs)

print(corr_1_2_5bcs)
print(corr_2_3_5bcs)
print(corr_3_1_5bcs)


# plot correlations
from scipy.stats import gaussian_kde

# Creating a 3x3 scatter plot grid
fig, axs = plt.subplots(3, 3, figsize=(15, 15))
fig.suptitle("RNA DNA ratio Correlation between reps", fontsize=16)

# Define datasets and titles for the plots
datasets = {
    "All Barcodes": rep_all_bcs,
    "10 Barcodes Threshold": rep10bcs,
    "5 Barcodes Threshold": rep5bcs
}

titles = ["Rep1 vs Rep2", "Rep2 vs Rep3", "Rep3 vs Rep1"]

# Populate scatter plots
for i, (label, dataset) in enumerate(datasets.items()):
    for j, title in enumerate(titles):
    
        x = dataset[j]["ratio_x"].values
        y = dataset[j]["ratio_y"].values

        # Calculate the point density
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)  # Density estimation

        # Scatter plot with density-based color
        scatter = axs[i, j].scatter(x, y, c=z, cmap='viridis', s=10, alpha=0.7)
    
        # Scatter plot
        #axs[i, j].scatter(dataset[j]["ratio_x"], dataset[j]["ratio_y"], alpha=0.7)
        # Labels and titles
        axs[i, j].set_title(f"{label} - {title}", fontsize=10)
        axs[i, j].set_xlabel("Ratio X", fontsize=8)
        axs[i, j].set_ylabel("Ratio Y", fontsize=8)
        axs[i, j].set_xlim(0, 20)  # Fix x-axis limit
        axs[i, j].set_ylim(0, 20)  # Fix y-axis limit
        axs[i, j].grid(True, alpha=0.3)

# Adjust layout and show the plot
plt.savefig(f'./additional/downsampling_rna_dna/scatters/downsampling{cells}_corr.png')