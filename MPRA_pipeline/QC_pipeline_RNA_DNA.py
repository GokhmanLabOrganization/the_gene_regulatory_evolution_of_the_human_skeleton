# NM 17-11-2024
# This script produces all the Relevant QC pipeline information and graphs 
# This script first requires all steps up to quantMPRA to be completed before running.

# Packages
import pandas as pd
from Bio import SeqIO
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import pickle



######################## Input ###############################################################
# 1. positive control oligos names
# 2. positive control oligos colors
# 3. negative control oligos names
# 4. negative control oligos colors
# 5. active oligos color



path = "/home/labs/davidgo/Collaboration/humanMPRA/additional/QC_pipeline/"

# Import positive control oligos list
pos_olg = pickle.load(open(path+"input/PosCtrl_osteoblast.pkl", "rb"))

# Import negative control oligos list 
neg_olg = pickle.load(open(path+"input/NegCtrl_not_active.pkl", "rb"))

# Import highlight oligos list
top_olg = pickle.load(open(path+"input/Hh_oligos.pkl", "rb"))

# Import df
activity_df = pd.read_csv(path+'input/comb_df_adjusted_fdr.csv')

# Color specifications
pos_col = "g"
neg_col = "r"
top_col = "y"

# Use CPM normalization?
cpm=True

# Use log scale in visualization?
logScale=True

#Additional filters?
min_DNA_reads = 5 


######################## Processing ###############################################################


# Add normalization
DNA_sum = activity_df["DNA_rep_comb"].sum()
RNA_sum = activity_df["RNA_rep_comb"].sum()

activity_df["DNA_rep_comb_cpm"] = 1000000*(activity_df["DNA_rep_comb"]+1)/DNA_sum
activity_df["RNA_rep_comb_cpm"] = 1000000*(activity_df["RNA_rep_comb"]+1)/RNA_sum

activity_df["DNA_rep_comb_cpm_log"] = np.log2(activity_df["DNA_rep_comb_cpm"])
activity_df["RNA_rep_comb_cpm_log"] = np.log2(activity_df["RNA_rep_comb_cpm"])


activity_df = activity_df[activity_df["DNA_rep_comb"] >= min_DNA_reads]

DNA_counts = "DNA_rep_comb"
RNA_counts = "RNA_rep_comb"

if cpm:
    DNA_counts=DNA_counts+"_cpm"
    RNA_counts=RNA_counts+"_cpm"

if logScale:
    DNA_counts=DNA_counts+"_log"
    RNA_counts=RNA_counts+"_log"




############################# RNA vs. DNA scatter plots ###########################

##### I. controls
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. log scale for both axes
# 2. background - colorized based on density
# 3. Positive controls colored with X color
# 4. Negative controls colored with Y color
# 5. Y=X line


#RNA DNA ratio
plt.clf()


# Prepare the data
x = activity_df[DNA_counts].values
y = activity_df[RNA_counts].values

# Create a 2D grid for the heatmap
x_grid, y_grid = np.meshgrid(
    np.linspace(x.min(), x.max(), 50),  # 500 points for high resolution
    np.linspace(y.min(), y.max(), 50)
)
grid_points = np.vstack([x_grid.ravel(), y_grid.ravel()])

# Evaluate the KDE on the grid
values = np.vstack([x, y])
kernel = gaussian_kde(values)
density = kernel(grid_points).reshape(x_grid.shape)

# Plot the heatmap
plt.figure(figsize=(8, 6))
plt.imshow(
    density,
    extent=[x.min(), x.max(), y.min(), y.max()],
    origin='lower',
    cmap='viridis',  # Change colormap if desired
    aspect='auto'
)

sns.scatterplot(data=activity_df[activity_df["oligo"].isin(pos_olg)], x=f'{DNA_counts}', y=f'{RNA_counts}', s=10,label = 'Positive control', color = pos_col)
sns.scatterplot(data=activity_df[activity_df["oligo"].isin(neg_olg)], x=f'{DNA_counts}', y=f'{RNA_counts}', s=10,label = 'Negative control', color = neg_col)

# Add Y = X line
x_limits = plt.xlim()
y_limits = plt.ylim()

min_val = max(x_limits[0], y_limits[0])
max_val = min(x_limits[1], y_limits[1])

plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, label='Y = X')


# Add labels or titles as needed
plt.title("RNA DNA graph with Controls")
plt.xlabel(DNA_counts)
plt.ylabel(RNA_counts)

plt.savefig(path+"output/RNA_DNA_ratio.png", dpi=1000)





##### II. Active oligos
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. Same as before
# 2. background - colorized based on density
# 3. color active oligos with Z color

plt.clf()

plt.figure(figsize=(8, 6))
plt.imshow(
    density,
    extent=[x.min(), x.max(), y.min(), y.max()],
    origin='lower',
    cmap='viridis',  # Change colormap if desired
    aspect='auto'
)

sns.scatterplot(data=activity_df[activity_df["oligo"].isin(top_olg)], x=f'{DNA_counts}', y=f'{RNA_counts}', s=10,label = 'Selected group', color = top_col)


# Add Y = X line
x_limits = plt.xlim()
y_limits = plt.ylim()

min_val = max(x_limits[0], y_limits[0])
max_val = min(x_limits[1], y_limits[1])

plt.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='--', linewidth=1, label='Y = X')


# Add labels or titles as needed
plt.title("RNA DNA graph with Controls")
plt.xlabel(DNA_counts)
plt.ylabel(RNA_counts)

plt.savefig(path+"output/selected_oligos_RNA_DNA_ratio.png", dpi=1000)


##### III. Per barcode
# df to use: ???

# important steps:
# 1. Same as before
# 2. instead of Oligos, print for each barcode before aggregation


############################# Between replicates correlations ##################################

# I. Active oligos
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. log scale for both axes
# 2. background - colorized based on density
# 3. Y=X line

# II. Selected group Active oligos
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. same as before
# 2. Mark a selected group of oligos


# III. Positive controls and Negative controls
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. same as before
# 2. Mark a selected group of oligos  
# 3. Positive controls colored with X color
# 4. Negative controls colored with Y color


# IV. RNA DNA ratio per replicate histogram
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. In the title - add the variance in each replicate (= the dynamic rannge of the replciate)


############################# Between batch correlations ######################################

# I. Positive controls Correlation
# df: humanMPRA\additional\analyze_controls\chondrocytes\controls_df.csv

# important steps:
# 1. Use the existing code, no further comments


############################# Activity plots ####################################################


# I. Active oligos
# df to use: LibraryAdaptor\output\activity_after_filter\comb_df_adjusted_fdr.csv

# important steps:
# 1. Use the existing code, no further comments



############################# Important statistics to print ####################################################
 
 
# Read counts summary statistics
print("Read counts summary statistics")

print(f"{'Type':<5} {'Median':<10} {'Mean':<10} {'Std Dev':<10} {'Median CPM':<15} {'Mean CPM':<15} {'Std CPM':<15}")

for data_type in ["DNA", "RNA"]:
    column = f"{data_type}_rep_comb"
    cpm_column = f"{data_type}_rep_comb_cpm"

    median = activity_df[column].median()
    mean = activity_df[column].mean()
    std = activity_df[column].std()
    median_cpm = activity_df[cpm_column].median()
    mean_cpm = activity_df[cpm_column].mean()
    std_cpm = activity_df[cpm_column].std()

    # Print results
    print(f"{data_type:<5} {median:<10.2f} {mean:<10.2f} {std:<10.2f} {median_cpm:<15.2f} {mean_cpm:<15.2f} {std_cpm:<15.2f}")
