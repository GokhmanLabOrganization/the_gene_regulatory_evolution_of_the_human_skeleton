# 14May2023 SF
# This script takes the output of quant_results.py and performs additional analyses and plots.

import pandas as pd
import sys
import os
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np
import re

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

# debug
library="L1"
adaptor="a2"
cells="chondrocytes"

minimal_barcodes = 5
output_figures_path = f'./{cells}/{library}{adaptor}/output/activity/outliers_analysis/'

comb_df_path = f'./{cells}/{library}{adaptor}/output/activity/comb_df.csv'
data_df = pd.read_csv(comb_df_path)
    
# add more columns to be used for additional calcs

# difference in log(RNA/DNA) between reps
data_df['rep2_rep1_delta_ratio_log'] = data_df['ratio_log_rep2'] - data_df['ratio_log_rep1']
data_df['rep3_rep1_delta_ratio_log'] = data_df['ratio_log_rep3'] - data_df['ratio_log_rep1']
data_df['rep2_rep3_delta_ratio_log'] = data_df['ratio_log_rep2'] - data_df['ratio_log_rep3']

# ratio between max RNA value and mean RNA value
data_df['rep1_max_RNA_vs_mean_RNA'] = data_df['max_UMI_RNA_rep1']/(data_df['RNA_pseudo_rep1']/data_df['possible_bcs'])
data_df['rep2_max_RNA_vs_mean_RNA'] = data_df['max_UMI_RNA_rep2']/(data_df['RNA_pseudo_rep2']/data_df['possible_bcs'])
data_df['rep3_max_RNA_vs_mean_RNA'] = data_df['max_UMI_RNA_rep3']/(data_df['RNA_pseudo_rep3']/data_df['possible_bcs'])

# standard deviation - of DNA and RNA counts per barcode
column_names=['DNA_' + cells + '_rep1','DNA_' + cells + '_rep2','DNA_' + cells + '_rep3','RNA_' + cells + '_rep1','RNA_' + cells + '_rep2','RNA_' + cells + '_rep3']
for column in column_names:
    for index, row in data_df.iterrows():
        data_df.loc[index,column + '_std'] = np.std([float(x) for x in list(data_df.loc[index,column].strip('[').strip(']').split(', '))])

# plots
sns.set(font_scale=0.35)

# create a filtered df, by minimal barcodes
data_df_filtered = data_df.copy()
data_df_filtered = data_df_filtered[(data_df_filtered["count_rep1"] >= minimal_barcodes) & (data_df_filtered["count_rep2"] >= minimal_barcodes) & (data_df_filtered["count_rep3"] >= minimal_barcodes)]

# fig 1: scatter plot between the difference in log(RNA/DNA) of two reps and outlier strength in one of them
for rep_pair in [["rep1", "rep2"], ["rep3", "rep2"], ["rep1", "rep3"]]:
    for rep in rep_pair:
        for activity in [0,1]:
            print(rep_pair)
            print(rep)
            print(activity)
            y_data_column = rep + '_max_RNA_vs_mean_RNA'
            x_data_column = rep_pair[1] + '_' + rep_pair[0] + "_delta_ratio_log"
            curr_data = data_df_filtered.copy()
            if activity == 1:
                curr_data = curr_data[curr_data["activity"] == "active"]
            values12 = np.vstack([curr_data[x_data_column], curr_data[y_data_column]])
            kernel12 = stats.gaussian_kde(values12)(values12)
            plt.clf()
            sns.scatterplot(curr_data, x=x_data_column, y=y_data_column, c=kernel12, s=5)
            plt.title(f'{cells}, {library}{adaptor}, max_RNA_vs_mean_RNA in {rep} vs delta {rep_pair[0]} and {rep_pair[1]} activity = {activity}')
            plt.ylabel('outlier_strength, ' + rep)
            plt.xlabel('delta RNA DNA ratio log, ' + rep_pair[0] + ' ' + rep_pair[1])
            plt.savefig(f'{output_figures_path}{cells}_{library}{adaptor}_max_RNA_vs_mean_RNA_outliers_{rep}_comparing_{rep_pair[0]}_{rep_pair[1]}_activity_{activity}.pdf')
            plt.savefig(f'{output_figures_path}{cells}_{library}{adaptor}_max_RNA_vs_mean_RNA_outliers_{rep}_comparing_{rep_pair[0]}_{rep_pair[1]}_activity_{activity}.png',dpi=300)

# average and median of log(RNA/DNA) differences
for rep_pair in [["rep1", "rep2"], ["rep3", "rep2"], ["rep1", "rep3"]]:
    for activity in [0, 1]:
        curr_data = data_df_filtered.copy()
        if activity == 1:
            curr_data = curr_data[curr_data["activity"] == "active"]
        average_diff = curr_data[rep_pair[1] + '_' + rep_pair[0] + "_delta_ratio_log"].mean()
        median_diff = curr_data[rep_pair[1] + '_' + rep_pair[0] + "_delta_ratio_log"].median()
        print('log(RNA/DNA) differences_' + rep_pair[0] + '_' + rep_pair[1] + '_' + str(activity) + '_average: ' + str(average_diff))
        print('log(RNA/DNA) differences_' + rep_pair[0] + '_' + rep_pair[1] + '_' + str(activity) + '_median: ' + str(median_diff))

# fig 2: scatter plot between the difference in log(RNA/DNA) of two reps and outlier strength in one of them - std of DNA or RNA
for rep_pair in [["rep1", "rep2"], ["rep3", "rep2"], ["rep1", "rep3"]]:
    for rep in rep_pair:
        for activity in [0,1]:
            for molecule in ['DNA', 'RNA']:
                print(rep_pair)
                print(rep)
                print(activity)
                print(molecule)
                y_data_column = molecule + '_' + cells + '_' + rep + '_std'
                x_data_column = rep_pair[1]+'_'+rep_pair[0]+"_delta_ratio_log"
                curr_data = data_df_filtered.copy()
                if activity == 1:
                    curr_data = curr_data[curr_data["activity"] == "active"]
                values12 = np.vstack([curr_data[x_data_column],curr_data[y_data_column]])
                kernel12 = stats.gaussian_kde(values12)(values12)
                plt.clf()
                sns.scatterplot(curr_data, x=x_data_column,y=y_data_column, c=kernel12, s=5)
                plt.title(f'{cells}, {library}{adaptor}, {molecule} _STD in {rep} vs delta {rep_pair[0]} and {rep_pair[1]} activity = {activity}')
                plt.ylabel(molecule + '_' + cells + '_' + rep + '_std')
                plt.xlabel('delta RNA DNA ratio log, ' + rep_pair[0] + ' ' + rep_pair[1])
                plt.savefig(f'{output_figures_path}{cells}_{library}{adaptor}_outliers_{molecule}_std_{rep}_comparing_{rep_pair[0]}_{rep_pair[1]}_activity_{activity}.pdf')
                plt.savefig(f'{output_figures_path}{cells}_{library}{adaptor}_outliers_{molecule}_std_{rep}_comparing_{rep_pair[0]}_{rep_pair[1]}_activity_{activity}.png',dpi=300)

# fig 3: histogram of DNA and RNA per barcode std
for molecule in ['DNA','RNA']:
    for rep in ["rep1", "rep2", "rep3"]:
        plt.clf()
        sns.displot(data=data_df_filtered, x=molecule + '_' + cells + '_' + rep + '_std', hue="activity", hue_order=["non_active", "active"],
                multiple="stack", palette=["lightgrey", "darkorange"])
        #sns.histplot(data=data_df_filtered, x=molecule + '_chondrocytes_' + rep + '_std', hue="activity")
        if molecule =='RNA':
            plt.xlim(0, 100)
        plt.savefig(f'{output_figures_path}{molecule}_std_per_barcode_{cells}_{library}{adaptor}_{rep}_color.pdf')
        plt.savefig(f'{output_figures_path}{molecule}_std_per_barcode_{cells}_{library}{adaptor}_{rep}_color.png',dpi=300)

# counts of average std
for molecule in ['DNA','RNA']:
    for rep in ["rep1", "rep2", "rep3"]:
        for activity in [0, 1]:
            curr_data = data_df_filtered.copy()
            if activity == 1:
                curr_data = curr_data[curr_data["activity"] == "active"]
            average_std = curr_data[molecule + '_' + cells + '_' + rep + '_std'].mean()
            median_std = curr_data[molecule + '_' + cells + '_' + rep + '_std'].median()
            print(molecule + '_' + cells + '_' + rep + '_activity' + str(activity) + '_std_average: ' + str(average_std))
            print(molecule + '_' + cells + '_' + rep + '_activity' + str(activity) + '_std_median: ' + str(median_std))

# fig 4: histogram of DNA and RNA count for all barcodes
for molecule in ['DNA','RNA']:
    for rep in ['rep1','rep2','rep3']:
        for activity in [0, 1]:



############################################### end of new analyses, the below was copied from the original script and is not used here ################

# boxplot
plt.clf()
sns.boxplot(data=data_df, x="alpha", y="control_type", order=['PosCtrl_neuron_active', 'PosCtrl_osteoblast_active', 'PosCtrl_chondrocyte_active', 'NegCtrl_active_not_diff', 'PosCtrl_diff', 'NegCtrl_not_active', 'scrambled_control', 'NegCtrl_non_SCREEN', 'No control'], palette=["gold", "gold", "gold", "gold", "gold", "silver", "silver", "silver", "lightcoral"])
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/alpha_per_control_type_boxplot_{cells}_{library}{adaptor}.pdf')

plt.clf()
sns.boxplot(data=data_df, x="alpha", y="control_type", showfliers = False, order=['PosCtrl_neuron_active', 'PosCtrl_osteoblast_active', 'PosCtrl_chondrocyte_active', 'NegCtrl_active_not_diff', 'PosCtrl_diff', 'NegCtrl_not_active', 'scrambled_control', 'NegCtrl_non_SCREEN', 'No control'], palette=["gold", "gold", "gold", "gold", "gold", "silver", "silver", "silver", "lightcoral"])
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/alpha_per_control_type_boxplot_{cells}_{library}{adaptor}_no_outliers.pdf')

# correlation between reps
corr_1_2_log = data_df["ratio_log_rep1"].corr(data_df["ratio_log_rep2"])
corr_2_3_log = data_df["ratio_log_rep2"].corr(data_df["ratio_log_rep3"])
corr_3_1_log = data_df["ratio_log_rep3"].corr(data_df["ratio_log_rep1"])

corr_1_2_10bcs_log = data_df.loc[(data_df["count_rep1"] >= 10)& (data_df["count_rep2"] >= 10), "ratio_log_rep1"].corr(data_df.loc[(data_df["count_rep1"] >= 10)& (data_df["count_rep2"] >= 10), "ratio_log_rep2"])
corr_2_3_10bcs_log = data_df.loc[(data_df["count_rep2"] >= 10)& (data_df["count_rep3"] >= 10), "ratio_log_rep2"].corr(data_df.loc[(data_df["count_rep2"] >= 10)& (data_df["count_rep3"] >= 10), "ratio_log_rep3"])
corr_3_1_10bcs_log = data_df.loc[(data_df["count_rep3"] >= 10)& (data_df["count_rep1"] >= 10), "ratio_log_rep3"].corr(data_df.loc[(data_df["count_rep3"] >= 10)& (data_df["count_rep1"] >= 10), "ratio_log_rep1"])

corr_1_2_5bcs_log = data_df.loc[(data_df["count_rep1"] >= 5)& (data_df["count_rep2"] >= 5), "ratio_log_rep1"].corr(data_df.loc[(data_df["count_rep1"] >= 5)& (data_df["count_rep2"] >= 5), "ratio_log_rep2"])
corr_2_3_5bcs_log = data_df.loc[(data_df["count_rep2"] >= 5)& (data_df["count_rep3"] >= 5), "ratio_log_rep2"].corr(data_df.loc[(data_df["count_rep2"] >= 5)& (data_df["count_rep3"] >= 5), "ratio_log_rep3"])
corr_3_1_5bcs_log = data_df.loc[(data_df["count_rep3"] >= 5)& (data_df["count_rep1"] >= 5), "ratio_log_rep3"].corr(data_df.loc[(data_df["count_rep3"] >= 5)& (data_df["count_rep1"] >= 5), "ratio_log_rep1"])

print("Correlation of Rep 1 and Rep2, log2: ", corr_1_2_log)
print("Correlation of Rep 2 and Rep3, log2: ", corr_2_3_log)
print("Correlation of Rep 3 and Rep1, log2: ", corr_3_1_log)
print("Correlation of Rep 1 and Rep2, log2 with >=10bcs: ", corr_1_2_10bcs_log)
print("Correlation of Rep 2 and Rep3, log2 with >=10bcs: ", corr_2_3_10bcs_log)
print("Correlation of Rep 3 and Rep1, log2 with >=10bcs: ", corr_3_1_10bcs_log)
print("Correlation of Rep 1 and Rep2, log2 with >=5bcs: ", corr_1_2_5bcs_log)
print("Correlation of Rep 2 and Rep3, log2 with >=5bcs: ", corr_2_3_5bcs_log)
print("Correlation of Rep 3 and Rep1, log2 with >=5bcs: ", corr_3_1_5bcs_log)


values12 = np.vstack([data_df[(data_df["count_rep1"] >= 5)& (data_df["count_rep2"] >= 5)]["ratio_log_rep1"], data_df[(data_df["count_rep1"] >= 5)& (data_df["count_rep2"] >= 5)]["ratio_log_rep2"]])
kernel12 = stats.gaussian_kde(values12)(values12)
plt.clf()
sns.scatterplot(data_df[(data_df["count_rep1"] >= 5)& (data_df["count_rep2"] >= 5)], x = "ratio_log_rep1", y = "ratio_log_rep2", c=kernel12, s=10)
plt.title(f'{cells}, {library}{adaptor}, correlation  rep1 and rep2 log2')
plt.xlabel('RNA/DNA rep1')
plt.ylabel('RNA/DNA rep2')
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_5bcs_rep12_log2_col_density.pdf')

values23 = np.vstack([data_df[(data_df["count_rep2"] >= 5)& (data_df["count_rep3"] >= 5)]["ratio_log_rep2"], data_df[(data_df["count_rep2"] >= 5)& (data_df["count_rep3"] >= 5)]["ratio_log_rep3"]])
kernel23 = stats.gaussian_kde(values23)(values23)
plt.clf()
sns.scatterplot(data_df[(data_df["count_rep2"] >= 5)& (data_df["count_rep3"] >= 5)], x = "ratio_log_rep2", y = "ratio_log_rep3", c=kernel23, s=10)
plt.title(f'{cells}, {library}{adaptor}, correlation  rep2 and rep3 log2')
plt.xlabel('RNA/DNA rep2')
plt.ylabel('RNA/DNA rep3')
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_5bcs_rep23_log2_col_density.pdf')

values31 = np.vstack([data_df[(data_df["count_rep3"] >= 5)& (data_df["count_rep1"] >= 5)]["ratio_log_rep3"], data_df[(data_df["count_rep3"] >= 5)& (data_df["count_rep1"] >= 5)]["ratio_log_rep1"]])
kernel31 = stats.gaussian_kde(values31)(values31)
plt.clf()
sns.scatterplot(data_df[(data_df["count_rep3"] >= 5)& (data_df["count_rep1"] >= 5)], x = "ratio_log_rep3", y = "ratio_log_rep1", c=kernel31, s=10)
plt.title(f'{cells}, {library}{adaptor}, correlation  rep3 and rep1 log2')
plt.xlabel('RNA/DNA rep3')
plt.ylabel('RNA/DNA rep1')
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_5bcs_rep31_log2_col_density.pdf')

# correlation for only active sequences
for rep in [["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]:
    print(rep)
    corr_active_all = data_df.loc[(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])
    corr_active_5 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])
    corr_active_10 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])

    print("all bcs:", corr_active_all)
    print(">=5bcs:", corr_active_5)
    print(">=10bcs:", corr_active_10)
    
    values = np.vstack([data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"ratio_log_{rep[0]}"], data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"ratio_log_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    plt.clf()
    sns.scatterplot(data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", c=kernel, s=10)
    plt.title(f'{cells}, {library}{adaptor}, correlation {rep[0]} and {rep[1]} rna/dna, log2 - active, >=5bcs')
    plt.xlabel(f'RNA/DNA, log2 {rep[0]}')
    plt.ylabel(f'RNA/DNA, log2 {rep[1]}')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_{rep[0]}_{rep[1]}_rna_dna.pdf')

#correlation plot of rna/dna - colored according to barcode counts
plt.clf()
fig, axes = plt.subplots(3, 3)
fig.suptitle('Correlation of RNA/DNA between replicates, colored by #bcs')
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="count_rep1", s=5,alpha = .9)
    sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="count_rep2", s=5,alpha = .9)
    sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="count_rep3", s=5,alpha = .9)
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_bcs_.pdf')

# correlation plot of rna/dna - colored according to dna counts
plt.clf()
fig, axes = plt.subplots(3, 3)
fig.suptitle('Correlation of RNA/DNA between replicates, colored by dna cpm')
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="DNA_pseudo_rep1", s=5,alpha = .9)
    sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="DNA_pseudo_rep2", s=5,alpha = .9)
    sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="DNA_pseudo_rep3", s=5,alpha = .9)
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_dna.pdf')

# correlation plot of rna/dna - colored according to rna counts
plt.clf()
fig, axes = plt.subplots(3, 3)
fig.suptitle('Correlation of RNA/DNA between replicates, colored by rna cpm')
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="RNA_pseudo_rep1", s=5,alpha = .9)
    sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="RNA_pseudo_rep2", s=5,alpha = .9)
    sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="RNA_pseudo_rep3", s=5,alpha = .9)
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_rna.pdf')

# correlation plot of rna/dna - colored according to maxUMI dna counts
plt.clf()
fig, axes = plt.subplots(3, 3)
fig.suptitle('Correlation of RNA/DNA between replicates, colored by maxUMIdna')
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_DNA_rep1", s=5,alpha = .9)
    sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_DNA_rep2", s=5,alpha = .9)
    sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_DNA_rep3", s=5,alpha = .9)
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_maxUMIdna.pdf')

# correlation plot of rna/dna - colored according to maxUMI rna counts
plt.clf()
fig, axes = plt.subplots(3, 3)
fig.suptitle('Correlation of RNA/DNA between replicates, colored by maxUMIrna')
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_RNA_rep1", s=5,alpha = .9)
    sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_RNA_rep2", s=5,alpha = .9)
    sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_RNA_rep3", s=5,alpha = .9)
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_maxUMIrna.pdf')

# correlation of rna cpm between replicates
print("correlation of RNA cpm between replicates")
for rep in [["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]:
    print(rep)
    corr_active_all = data_df.loc[(data_df["activity"] == "active"), f"RNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df["activity"] == "active"), f"RNA_pseudo_{rep[1]}"])
    corr_active_5 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[1]}"])
    corr_active_10 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[1]}"])

    print("all bcs:", corr_active_all)
    print(">=5bcs:", corr_active_5)
    print(">=10bcs:", corr_active_10)
    
    
    values = np.vstack([data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"RNA_pseudo_{rep[0]}"], data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"RNA_pseudo_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    plt.clf()
    sns.scatterplot(data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"RNA_pseudo_{rep[0]}", y = f"RNA_pseudo_{rep[1]}", c=kernel, s=10)
    plt.title(f'{cells}, {library}{adaptor}, correlation {rep[0]} and {rep[1]} RNA cpm - active, >=5bcs')
    plt.xlabel(f'RNA cpm {rep[0]}')
    plt.ylabel(f'RNA cpm {rep[1]}')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_{rep[0]}_{rep[1]}_RNA.pdf')

# correlation of dna cpm between replicates
print("correlation of DNA cpm between replicates")
for rep in [["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]:
    print(rep)
    corr_active_all = data_df.loc[(data_df["activity"] == "active"), f"DNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df["activity"] == "active"), f"DNA_pseudo_{rep[1]}"])
    corr_active_5 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[1]}"])
    corr_active_10 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[1]}"])

    print("all bcs:", corr_active_all)
    print(">=5bcs:", corr_active_5)
    print(">=10bcs:", corr_active_10)
    
    
    values = np.vstack([data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"DNA_pseudo_{rep[0]}"], data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"DNA_pseudo_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    plt.clf()
    sns.scatterplot(data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"DNA_pseudo_{rep[0]}", y = f"DNA_pseudo_{rep[1]}", c=kernel, s=10)
    plt.title(f'{cells}, {library}{adaptor}, correlation {rep[0]} and {rep[1]} DNA cpm - active, >=5bcs')
    plt.xlabel(f'DNA cpm {rep[0]}')
    plt.ylabel(f'DNA cpm {rep[1]}')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_active_5bcs_{rep[0]}_{rep[1]}_DNA.pdf')

sns.reset_defaults()

# histogram of RNA/DNA colored according to activity and stacked
for rep in ["rep1", "rep2", "rep3"]:
    plt.clf()
    sns.displot(data=data_df, x=f"ratio_log_{rep}", hue = "activity", hue_order = ["non_active", "active"], multiple="stack", palette=["lightgrey", "darkorange"])
    plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/RNA_DNA_distribution_stacked_{cells}_{library}{adaptor}_{rep}_color.pdf')

# distribution of alphas
plt.clf()
sns.histplot(data=data_df, x="alpha", hue = "activity")
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/alpha_distribution_{cells}_{library}{adaptor}.pdf')

plt.clf()
sns.histplot(data=data_df, x="alpha", hue = "activity")
plt.xlim(0, 20)
plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/alpha_distribution_{cells}_{library}{adaptor}_xlim_20.pdf')

# distribution of scrambled
for rep in ["rep1", "rep2", "rep3"]:
    plt.clf()
    data_scrambled = data_df[data_df["control_type"]== "scrambled_control"]
    sns.histplot(data=data_scrambled, x=f"ratio_log_{rep}")
    plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/scrambled_distribution_{cells}_{library}{adaptor}_{rep}.pdf')
    
# get oligo coverage for each replicate
print("oligo coverage")
for rep in ["rep1", "rep2", "rep3"]:
    for bc in [1, 5, 10, 20]:
        count = (data_df[f"count_{rep}"] >= bc).sum()
        ratio = count/79812
        print(f"{rep}: count of oligos with >= {bc} barcodes: {count}")
        print(f"{rep}: ratio of oligos with >= {bc} barcodes: {ratio}")
