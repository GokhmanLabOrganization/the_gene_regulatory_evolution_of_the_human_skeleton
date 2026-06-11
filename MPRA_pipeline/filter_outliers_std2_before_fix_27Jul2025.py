# the aim of this script is to apply an outlier filter (2 stds) on barcode_counts_UMI.txt and return a new file in the format of barcode_counts_UMI.txt that can then be used as input for the rest of the pipeline.
# in addition the script should also check the effect of the filter (e.g. how many barcodes are removed, correlation etc)

import pandas as pd
import numpy as np
import sys
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import os
import seaborn as sns
from collections import Counter
from statistics import mean
import pickle


print(matplotlib.__version__)

###################### Changes introduced Simon 28May2025:
# 1. Fix the bug that caused oligo-bc pairs being outliers in all 3 replicates in parallel to not being filtered out
# 2. A new output to log the filtering - a table with filtered out oligo-bc pairs and the number of replicates they are outliers in
# Most changes were done in the core function std_n_outliers
# Previously: for each replicate the oligo_bc pairs with zero RNA and zero DNA were filtered before entering the function.
# Therefore each replicate df had a different number of lines, less or equal to the original table with all data.
# Each replicate, after filtering within the function, was smaller than the original table, and null values were replaced with 0 after the function.
# The bug was when trying to fill NA data with zeros, instead the script was copying the raw data from the raw table for lines having outliers in all 3 replicates, keeping the outlier data.
# Now: the entire replicate table enters the function and leaves the function with same number of rows. Zeros replace outliers within the function.
######################

#todo
#add counts for other libraries x
#delete other filter adn barcode thresholds x
#save exploded df x
#statistics of filter: histogram of how many bcs per oligo are removed, oligo coverage, correlation between reps + average correlation between reps (for different dna count thresholds and activity thresholds)
#scatterpliot of correlation

def std_n_outliers(rna, dna, oligo_bc, n):
    mean = np.mean(rna)
    std = np.std(rna)
    lower_bound = mean - n * std
    upper_bound = mean + n * std
    indices = [i for i, v in enumerate(rna) if lower_bound < v < upper_bound]
    rna = [rna[i] for i in indices]
    dna = [dna[i] for i in indices]
    oligo_bc = [oligo_bc[i] for i in indices]
    return rna, dna, oligo_bc

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

def create_ratio_df():


    # access replicate UMI counts files
    rep1 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep1.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])
    rep2 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep2.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])
    rep3 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep3.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])

    UMI_all_reps = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI.txt', sep='\t')
    UMI_all_reps.set_index('oligo_bc', inplace=True)

    #Obtaining the read counts with soft coding (replacing the numbers above), Nadav 29-0-2024
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
    

    # get number of reads 
    DNA_reads = [rep1DNA, rep2DNA, rep3DNA]
    RNA_reads = [rep1RNA, rep2RNA, rep3RNA]
    read_dict = {"DNA_reads": DNA_reads, "RNA_reads": RNA_reads}
    print(DNA_reads)
    print(RNA_reads)

    # group oligos and sum RNA, DNA counts. Also counts how many barcodes we have per oligo. Calculate ratio and perform log2 on ratio
    replicates = []
    exploded_df = []
    
    for i, rep in enumerate([rep1, rep2, rep3]):
        rep[f'oligo_bc_filtered_std2'] = None
        for molecule in ["RNA", "DNA"]:
            rep[f'{molecule}_filtered_std2'] = None
        rep = rep[(rep['DNA_no_filter'] > 0) | (rep['RNA_no_filter'] > 0)]
        rep[['oligo', 'bc']] = rep.oligo_bc_no_filter.str.rsplit("_", n=1, expand=True) #30.09.2024 NM - changed  rsplit("_", 1, expand=True) to rsplit("_", n=1, expand=True)
        rep = rep.set_index('oligo')
        rep_list = rep.groupby('oligo').aggregate(lambda a: a.tolist())
        rep_list = rep_list.reset_index()
        for index, row in rep_list.iterrows():
            rep_list.at[index, 'RNA_filtered_std2'], rep_list.at[index, 'DNA_filtered_std2'], rep_list.at[index, 'oligo_bc_filtered_std2'] = std_n_outliers(
                row['RNA_no_filter'], row['DNA_no_filter'], row['oligo_bc_no_filter'], 2)
        for outlier_filter in ["no_filter","filtered_std2"]:
            to_explode = rep_list[[col for col in rep_list.columns if f"{outlier_filter}" in col]].copy()
            to_explode = to_explode[to_explode[f'oligo_bc_{outlier_filter}'].astype(bool)]
            exploded = to_explode.explode([f"RNA_{outlier_filter}", f"DNA_{outlier_filter}", f'oligo_bc_{outlier_filter}'], ignore_index=True)
            exploded_df.append(exploded)
            rep_list = rep_list.drop([f'oligo_bc_{outlier_filter}'], axis=1)
        rep_list = rep_list.drop(['bc'], axis=1)
        for outlier_filter in ["no_filter", "filtered_std2"]:
            rep_list[f'count_{outlier_filter}'] = rep_list[f'RNA_{outlier_filter}'].apply(len)
            mask = (rep_list[f'count_{outlier_filter}'] > 0) # filtering might lead to some cases where a barcode ends up with 0 barcodes - want to filter out these oligos
            rep_valid = rep_list[mask] # This variable is not used anywhere. Perhaps this causes the bug?  # Nadav 24.06.2024
            for molecule in ["RNA", "DNA"]:
                rep_list.loc[mask, f'{molecule}_{outlier_filter}_sum'] = rep_list.loc[mask, f'{molecule}_{outlier_filter}'].apply(sum)
                rep_list.loc[mask, f"{molecule}_{outlier_filter}_pseudo"] = rep_list.loc[mask, f"{molecule}_{outlier_filter}_sum"] + 1
                rep_list.loc[mask, f"{molecule}_{outlier_filter}_cpm"] = (rep_list.loc[mask, f"{molecule}_{outlier_filter}_pseudo"] * 1000000)/read_dict[f"{molecule}_reads"][0]  # todo: this has to be fixed - change 0 to i

            rep_list.loc[mask, f"ratio_{outlier_filter}"] = rep_list.loc[mask, f"RNA_{outlier_filter}_cpm"] / rep_list.loc[mask, f"DNA_{outlier_filter}_cpm"]
            rep_list.loc[mask, f"ratio_log_{outlier_filter}"] = np.log2(rep_list.loc[mask, f"ratio_{outlier_filter}"])
        replicates.append(rep_list)

    rep1_counts = replicates[0]
    rep2_counts = replicates[1]
    rep3_counts = replicates[2]

    # merge the different outputs
    ratio_df = rep1_counts.merge(rep2_counts,on='oligo',how="outer",suffixes=("_rep1",None)).merge(rep3_counts,on='oligo',how="outer",suffixes=("_rep2","_rep3"))
    # save df
    ratio_df.to_csv(f'./{cells}/{library}{adaptor}/output/filter/ratio_wo_outliers_std2.csv', header=True, index=False)
    
    # save the exploded df (2std) - so we get a file in the same format as barcode_counts_UMI.txt - which than can be used to create the inout for MPRAnalyze and the other steps of the analysis (get one file for all the replicates)
    

    rep1_exploded_filtered = exploded_df[1]
    rep2_exploded_filtered = exploded_df[3]
    rep3_exploded_filtered = exploded_df[5]
    exploded_df_filter = rep1_exploded_filtered.merge(rep2_exploded_filtered,on='oligo_bc_filtered_std2',how="outer",suffixes=("_rep1",None)).merge(rep3_exploded_filtered,on='oligo_bc_filtered_std2',how="outer",suffixes=("_rep2","_rep3"))
    # print(list(exploded_df_filter.columns))
    exploded_df_filter.fillna(0.0, inplace=True)
    exploded_df_filter = exploded_df_filter[["oligo_bc_filtered_std2", "RNA_filtered_std2_rep1", "RNA_filtered_std2_rep2","RNA_filtered_std2_rep3","DNA_filtered_std2_rep1","DNA_filtered_std2_rep2","DNA_filtered_std2_rep3"]]
    exploded_df_filter.columns = ["oligo_bc",f"RNA_{cells}_rep1",f"RNA_{cells}_rep2",f"RNA_{cells}_rep3",f"DNA_{cells}_rep1",f"DNA_{cells}_rep2",f"DNA_{cells}_rep3"]
    exploded_df_filter.set_index('oligo_bc', inplace=True)
    UMI_comb_filter = exploded_df_filter.combine_first(UMI_all_reps)
    UMI_comb_filter.reset_index(inplace=True)
    UMI_comb_filter.to_csv(f'./{cells}/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', index=False, sep='\t')

    rep1_exploded_no_filter = exploded_df[0]
    rep2_exploded_no_filter = exploded_df[2]
    rep3_exploded_no_filter = exploded_df[4]
    
    exploded_df_no_filter = rep1_exploded_no_filter.merge(rep2_exploded_no_filter,on='oligo_bc_no_filter',how="outer",suffixes=("_rep1",None)).merge(rep3_exploded_no_filter,on='oligo_bc_no_filter',how="outer",suffixes=("_rep2","_rep3"))
    exploded_df_no_filter.fillna(0.0, inplace=True)
    exploded_df_no_filter = exploded_df_no_filter[["oligo_bc_no_filter", "RNA_no_filter_rep1", "RNA_no_filter_rep2","RNA_no_filter_rep3","DNA_no_filter_rep1","DNA_no_filter_rep2","DNA_no_filter_rep3"]]
    exploded_df_no_filter.columns = ["oligo_bc",f"RNA_{cells}_rep1",f"RNA_{cells}_rep2",f"RNA_{cells}_rep3",f"DNA_{cells}_rep1",f"DNA_{cells}_rep2",f"DNA_{cells}_rep3"]
    exploded_df_no_filter.set_index('oligo_bc', inplace=True)
    UMI_comb_no_filter = exploded_df_no_filter.combine_first(UMI_all_reps)
    UMI_comb_no_filter.reset_index(inplace=True)
    UMI_comb_no_filter.to_csv(f'./{cells}/{library}{adaptor}/output/UMI/UMI_exploded_no_filter.txt', index=False, sep='\t')

    return ratio_df

if not os.path.exists(f'./{cells}/{library}{adaptor}/output/filter/ratio_wo_outliers_std2.csv'):
    ratio_df = create_ratio_df()
else:
    ratio_df = pd.read_csv(f'./{cells}/{library}{adaptor}/output/filter/ratio_wo_outliers_std2.csv')
    
    
# overview over how many we loose per oligo for the different methods and the difference in RNA/DNA before and after filter (for scrambked and for all)
   
plt.clf()
fig, axes = plt.subplots(6, 2, sharex="col")
for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    for m, column in enumerate(["count","ratio_log"]):
            print(rep,column)
            if column=="count":
                num_bins= int((ratio_df[f"count_no_filter_{rep}"]-ratio_df[f"count_filtered_std2_{rep}"]).max() + 1)
                mini= 0
                maxi= (ratio_df[f"count_no_filter_{rep}"]-ratio_df[f"count_filtered_std2_{rep}"]).max()
            else:
                num_bins=60
                mini=(ratio_df[f"ratio_log_no_filter_{rep}"]-ratio_df[f"ratio_log_filtered_std2_{rep}"]).min()
                maxi=(ratio_df[f"ratio_log_no_filter_{rep}"]-ratio_df[f"ratio_log_filtered_std2_{rep}"]).max()
            axes[n,m].hist(ratio_df[f"{column}_no_filter_{rep}"]-ratio_df[f"{column}_filtered_std2_{rep}"], bins=num_bins, range=(mini,maxi)) #, bins=25, range=(0,24))
            scrambled_df = ratio_df[ratio_df['oligo'].str.contains("scrambled")].copy()
            axes[n+3,m].hist(scrambled_df[f"{column}_no_filter_{rep}"]-scrambled_df[f"{column}_filtered_std2_{rep}"], bins=num_bins, range=(mini,maxi)) #, bins=25, range=(0,24))
axes[-1,0].set(xlabel="count - count 2std")
axes[-1,1].set(xlabel="ratio - ratio 2std")
        # print(Counter(ratio_df[f"count_no_filter_{rep}"]-ratio_df[f"count_{outlier_filter}_{rep}"]))
plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/histogram_number_ratio_filter.pdf')


# correlation between the reps for the different methods and different dna count thresholds and activity thresholds

for outlier_filter in ["no_filter","filtered_std2"]:
    plt.clf()
    fig, axes = plt.subplots(4,3, sharey=True, sharex=True)
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    x_label_list = ['rep1', 'rep2', 'rep3']
    y_label_list = ['rep1', 'rep2', 'rep3']
    for n, activity_level in enumerate(["none", 0, 1]):
            for m, threshold in enumerate([0, 5, 10, 20]):
                for rep in ['rep1', 'rep2', 'rep3']:
                    ratio_df[f'ratio_{outlier_filter}_{rep}_DNA_{threshold}'] = ratio_df[f"ratio_log_{outlier_filter}_{rep}"].where(ratio_df[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold, pd.NA) # only log ratio where DNA count of rep is higher than threshold
                if activity_level=="none":
                    ratio_df_copy = ratio_df.copy()
                else:
                    ratio_df_copy = ratio_df[(ratio_df[f'ratio_{outlier_filter}_rep1_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep2_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep3_DNA_{threshold}'] > activity_level)].copy()
                corr = ratio_df_copy[[f'ratio_{outlier_filter}_rep1_DNA_{threshold}',f'ratio_{outlier_filter}_rep2_DNA_{threshold}',f'ratio_{outlier_filter}_rep3_DNA_{threshold}']].corr()
                if activity_level =="none" and threshold ==5:
                    avg_corr = mean([corr.iloc[1,0], corr.iloc[2,0], corr.iloc[2,1]])
                    print(f"avergae correlation for {outlier_filter}, activity threshold: {activity_level}, dna counts threshold: {threshold} - {avg_corr}")
                sns.heatmap(ax=axes[m,n], data=corr, cmap="Reds", vmin=0, vmax=1, annot=True, annot_kws={"fontsize":6}, fmt=".2f", square=True, xticklabels=x_label_list, yticklabels=y_label_list, cbar_ax=cbar_ax)
                axes[m, n].tick_params(axis='x', labelrotation = 90)
                axes[-1, n].set(xlabel=f'{activity_level}')
                axes[m, 0].set(ylabel=f'{threshold}')
    # fig.supxlabel('activity threshold')
    # fig.supylabel('DNA threshold')
    plt.tight_layout()      
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/correlation_between_reps_{outlier_filter}_and_DNA_filter_activity_filter.pdf')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/correlation_between_reps_{outlier_filter}_and_DNA_filter_activity_filter.png', dpi=1000)

# without top 1 percentile oligos with highest DNA counts 
filter_cols = ["DNA_filtered_std2_sum_rep1","DNA_filtered_std2_sum_rep2","DNA_filtered_std2_sum_rep3","RNA_filtered_std2_sum_rep1","RNA_filtered_std2_sum_rep2","RNA_filtered_std2_sum_rep3"] #for number of DNA counts
#filter_cols = ["count_filtered_std2_rep1","count_filtered_std2_rep2","count_filtered_std2_rep3"] # for number of DNA barcodes
#filter_cols = ["count_filtered_std2_rep1","count_filtered_std2_rep2","count_filtered_std2_rep3","DNA_filtered_std2_sum_rep1","DNA_filtered_std2_sum_rep2","DNA_filtered_std2_sum_rep3","RNA_filtered_std2_sum_rep1","RNA_filtered_std2_sum_rep2","RNA_filtered_std2_sum_rep3"] # for number of DNA barcodes


percentiles = {col: ratio_df[col].quantile(0.95) for col in filter_cols}
print("percentiles:",percentiles)
ratio_df_99_per = ratio_df[(ratio_df[filter_cols] < pd.Series(percentiles)).all(axis=1)]
print(len(ratio_df))
print(len(ratio_df_99_per))

for outlier_filter in ["no_filter","filtered_std2"]:
    plt.clf()
    fig, axes = plt.subplots(4,3, sharey=True, sharex=True)
    cbar_ax = fig.add_axes([.91, .3, .03, .4])
    x_label_list = ['rep1', 'rep2', 'rep3']
    y_label_list = ['rep1', 'rep2', 'rep3']
    for n, activity_level in enumerate(["none", 0, 1]):
            for m, threshold in enumerate([0, 5, 10, 20]):
                for rep in ['rep1', 'rep2', 'rep3']:
                    ratio_df_99_per[f'ratio_{outlier_filter}_{rep}_DNA_{threshold}'] = ratio_df_99_per[f"ratio_log_{outlier_filter}_{rep}"].where(ratio_df_99_per[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold, pd.NA) # only log ratio where DNA count of rep is higher than threshold
                if activity_level=="none":
                    ratio_df_99_per_copy = ratio_df_99_per.copy()
                else:
                    ratio_df_99_per_copy = ratio_df_99_per[(ratio_df_99_per[f'ratio_{outlier_filter}_rep1_DNA_{threshold}'] > activity_level) |(ratio_df_99_per[f'ratio_{outlier_filter}_rep2_DNA_{threshold}'] > activity_level) |(ratio_df_99_per[f'ratio_{outlier_filter}_rep3_DNA_{threshold}'] > activity_level)].copy()
                corr = ratio_df_99_per_copy[[f'ratio_{outlier_filter}_rep1_DNA_{threshold}',f'ratio_{outlier_filter}_rep2_DNA_{threshold}',f'ratio_{outlier_filter}_rep3_DNA_{threshold}']].corr()
                if activity_level =="none" and threshold ==5:
                    avg_corr = mean([corr.iloc[1,0], corr.iloc[2,0], corr.iloc[2,1]])
                    print(f"avergae correlation for {outlier_filter}, activity threshold: {activity_level}, dna counts threshold: {threshold} - {avg_corr}")
                sns.heatmap(ax=axes[m,n], data=corr, cmap="Reds", vmin=0, vmax=1, annot=True, annot_kws={"fontsize":6}, fmt=".2f", square=True, xticklabels=x_label_list, yticklabels=y_label_list, cbar_ax=cbar_ax)
                axes[m, n].tick_params(axis='x', labelrotation = 90)
                axes[-1, n].set(xlabel=f'{activity_level}')
                axes[m, 0].set(ylabel=f'{threshold}')
    # fig.supxlabel('activity threshold')
    # fig.supylabel('DNA threshold')
    plt.tight_layout()      
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/correlation_between_reps_{outlier_filter}_and_DNA_filter_activity_filter_99_percentile.pdf')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/correlation_between_reps_{outlier_filter}_and_DNA_filter_activity_filter_99_percentile.png', dpi=1000)

# scatter plots of correlation between rep1 and rep2 for the different methods and dna count thresholds and activity thresholds
for outlier_filter in ["no_filter","filtered_std2"]:
    plt.clf()
    fig, axes = plt.subplots(4,3, sharey=True, sharex=True)
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                ratio_df_99_per[f'ratio_{outlier_filter}_{rep}_DNA_{threshold}'] = ratio_df_99_per[f"ratio_log_{outlier_filter}_{rep}"].where(ratio_df_99_per[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold, pd.NA)
            if activity_level=="none":
                ratio_df_99_per = ratio_df_99_per.copy()
            else:
                ratio_df_99_per = ratio_df_99_per[(ratio_df_99_per[f'ratio_{outlier_filter}_rep1_DNA_{threshold}'] > activity_level) |(ratio_df_99_per[f'ratio_{outlier_filter}_rep2_DNA_{threshold}'] > activity_level) |(ratio_df_99_per[f'ratio_{outlier_filter}_rep3_DNA_{threshold}'] > activity_level)].copy()
            sns.scatterplot(ax=axes[m, n], data = ratio_df_99_per, x = f'ratio_{outlier_filter}_rep1_DNA_{threshold}', y = f'ratio_{outlier_filter}_rep2_DNA_{threshold}', s=5)

            axes[-1, n].set(xlabel=f'{activity_level}')
            axes[m, 0].set(ylabel=f'{threshold}')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/scatter_correlation_between_reps_{outlier_filter}_and_DNA_filter_activityfilter_99_percentile.pdf')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/scatter_correlation_between_reps_{outlier_filter}_and_DNA_filter_activityfilter_99_percentile.png', dpi=1000)





# delta corr

plt.clf()
fig, axes = plt.subplots(4,3, sharey=True, sharex=True)
cbar_ax = fig.add_axes([.91, .3, .03, .4])
x_label_list = ['rep1', 'rep2', 'rep3']
y_label_list = ['rep1', 'rep2', 'rep3']
for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                ratio_df[f'ratio_no_filter_{rep}_DNA_{threshold}'] = ratio_df[f"ratio_log_no_filter_{rep}"].where(ratio_df[f"DNA_no_filter_sum_{rep}"] >= threshold, pd.NA) # only log ratio where DNA count of rep is higher than threshold
                ratio_df[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = ratio_df[f"ratio_log_filtered_std2_{rep}"].where(ratio_df[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA) # only log ratio where DNA count of rep is higher than threshold
            if activity_level=="none":
                ratio_df_copy = ratio_df.copy()
            else:
                ratio_df_copy_no_filter = ratio_df[(ratio_df[f'ratio_no_filter_rep1_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_no_filter_rep2_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_no_filter_rep3_DNA_{threshold}'] > activity_level)].copy()
                ratio_df_copy_filter = ratio_df[(ratio_df[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
            corr_filter = ratio_df_copy[[f'ratio_filtered_std2_rep1_DNA_{threshold}',f'ratio_filtered_std2_rep2_DNA_{threshold}',f'ratio_filtered_std2_rep3_DNA_{threshold}']].corr()
            corr_no_filter = ratio_df_copy[[f'ratio_no_filter_rep1_DNA_{threshold}',f'ratio_no_filter_rep2_DNA_{threshold}',f'ratio_no_filter_rep3_DNA_{threshold}']].corr()
            delta_corr = corr_filter - corr_no_filter.values
            # print(delta_corr)
            sns.heatmap(ax=axes[m,n], data=delta_corr, cmap="Reds", annot=True, annot_kws={"fontsize":6}, fmt=".2f", square=True, xticklabels=x_label_list, yticklabels=y_label_list, cbar_ax=cbar_ax)
            axes[m, n].tick_params(axis='x', labelrotation = 90)
            axes[-1, n].set(xlabel=f'{activity_level}')
            axes[m, 0].set(ylabel=f'{threshold}')
# fig.supxlabel('activity threshold')
# fig.supylabel('DNA threshold')
plt.tight_layout()
plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/delta_correlation_between_reps_and_DNA_filter_activity_filter.pdf')


# scatter plots of correlation between rep1 and rep2 for the different methods and dna count thresholds and activity thresholds
for outlier_filter in ["no_filter","filtered_std2"]:
    plt.clf()
    fig, axes = plt.subplots(4,3, sharey=True, sharex=True)
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                ratio_df[f'ratio_{outlier_filter}_{rep}_DNA_{threshold}'] = ratio_df[f"ratio_log_{outlier_filter}_{rep}"].where(ratio_df[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold, pd.NA)
            if activity_level=="none":
                ratio_df_copy = ratio_df.copy()
            else:
                ratio_df_copy = ratio_df[(ratio_df[f'ratio_{outlier_filter}_rep1_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep2_DNA_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep3_DNA_{threshold}'] > activity_level)].copy()
            sns.scatterplot(ax=axes[m, n], data = ratio_df_copy, x = f'ratio_{outlier_filter}_rep1_DNA_{threshold}', y = f'ratio_{outlier_filter}_rep2_DNA_{threshold}', s=5)

            axes[-1, n].set(xlabel=f'{activity_level}')
            axes[m, 0].set(ylabel=f'{threshold}')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/scatter_correlation_between_reps_{outlier_filter}_and_DNA_filter_activityfilter.pdf')
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/scatter_correlation_between_reps_{outlier_filter}_and_DNA_filter_activityfilter.png', dpi=1000)


# # calculate oligo coverage
## contains a bug!!! does not consider activity threshold
print("oligo coverage")
for outlier_filter in ["no_filter","filtered_std2", "delta"]:
    print(outlier_filter)
    for n, activity_level in enumerate(["none", 0, 1]):
        print(activity_level)
        for m, threshold in enumerate([0, 5, 10, 20]):
            print(threshold)
            for rep in ["rep1", "rep2", "rep3"]:
                if outlier_filter == "delta":
                    count_no_filter = (ratio_df[f"DNA_no_filter_sum_{rep}"] >= threshold).sum()
                    count_filter = (ratio_df[f"DNA_filtered_std2_sum_{rep}"] >= threshold).sum()
                    ratio_delta = round((count_no_filter - count_filter)/79812, 2)
                    print(rep, ratio_delta)
                else:
                    count = (ratio_df[f"DNA_{outlier_filter}_sum_{rep}"] >= threshold).sum()
                    ratio = round(count/79812, 2)
                    print(rep, ratio)
            
            



######################Group Analysis######################

#import the group names dictionary which contains the relevant group names for the specific library
oligo_groups= pickle.load( open( f'./{cells}/{library}{adaptor}/input/oligo_groups.pkl', "rb" ) )

#RNA DNA ratio

newpath = f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA_ratio/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
for key,val in oligo_groups.items():
    fig, axes = plt.subplots(4, 3, sharey=True, sharex=True)
    
    # Filter the DataFrame into two parts: one for the selected rows and one for the rest
    sub_group = ratio_df[ratio_df['oligo'].str.contains(key, na=False)].copy()
    non_group = ratio_df[~ratio_df['oligo'].str.contains(key, na=False)].copy()
    if key == "a3" or key == "a4":
        sub_group = ratio_df[
            ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        non_group = ratio_df[
            ~ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
    # Loop through the activity levels and thresholds to create subplots
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                non_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = non_group[f"ratio_log_filtered_std2_{rep}"].where(non_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA)
                sub_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = sub_group[f"ratio_log_filtered_std2_{rep}"].where(sub_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA)
                
                if activity_level=="none":
                    non_group_copy = non_group.copy()
                    sub_group_copy = sub_group.copy()
                else:
                    non_group_copy = non_group[(non_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(non_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(non_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    sub_group_copy = sub_group[(sub_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(sub_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(sub_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    
            sns.scatterplot(ax=axes[m, n], data=non_group_copy, x=f'ratio_filtered_std2_rep1_DNA_{threshold}', y=f'ratio_filtered_std2_rep2_DNA_{threshold}', s=5, color="gray", alpha=0.5,label=f'N={len(non_group_copy)}')
            sns.scatterplot(ax=axes[m, n], data=sub_group_copy, x=f'ratio_filtered_std2_rep1_DNA_{threshold}', y=f'ratio_filtered_std2_rep2_DNA_{threshold}', s=5, color="red",label=f'N={len(sub_group_copy)}')
            axes[m, n].legend(fontsize="2")
            axes[-1, n].set_xlabel(f'{activity_level}',fontsize=10)
            axes[m, 0].set_ylabel(f'{threshold}',fontsize=10)
            
            #Add correaltion
            sub_group_corr = sub_group_copy[[f'ratio_filtered_std2_rep1_DNA_{threshold}', f'ratio_filtered_std2_rep2_DNA_{threshold}']].corr().iloc[0, 1]
            if pd.notna(sub_group_corr):
                axes[m, n].text(0.05, 0.85, f'r={sub_group_corr:.2f}', transform=axes[m, n].transAxes, ha='left', va='top', fontsize=3, color="red")
            else:
                axes[m, n].text(0.05, 0.85, 'r=N/A', transform=axes[m, n].transAxes, ha='left', va='top', fontsize=3, color="red")


    fig.suptitle(f'Corr between reps 1 and 2 RNA DNA ratio\n{val}', fontsize=12)
    fig.supxlabel('Activity threshold',fontsize=7)
    fig.supylabel('Min DNA count',fontsize=7)
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA_ratio/RNA_DNA_ratio_corr_{val}.png', dpi=1000)

    plt.clf()

for key, val in oligo_groups.items():
    fig, axes = plt.subplots(4, 3, sharey=True, sharex=True)

    # Filter the DataFrame into two parts: one for the selected rows and one for the rest
    sub_group = ratio_df[ratio_df['oligo'].str.contains(key, na=False)].copy()
    non_group = ratio_df[~ratio_df['oligo'].str.contains(key, na=False)].copy()
    if key == "a3" or key == "a4":
        sub_group = ratio_df[
            ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        non_group = ratio_df[
            ~ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        # Loop through the activity levels and thresholds to create subplots
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                non_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = non_group[
                    f"ratio_log_filtered_std2_{rep}"].where(non_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold,
                                                            pd.NA)
                sub_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = sub_group[
                    f"ratio_log_filtered_std2_{rep}"].where(sub_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold,
                                                            pd.NA)

                if activity_level == "none":
                    non_group_copy = non_group.copy()
                    sub_group_copy = sub_group.copy()
                else:
                    non_group_copy = non_group[
                        (non_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) | (
                                    non_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) | (
                                    non_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    sub_group_copy = sub_group[
                        (sub_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) | (
                                    sub_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) | (
                                    sub_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()

            sns.histplot(ax=axes[m, n], data=non_group_copy, x=f'ratio_filtered_std2_rep1_DNA_{threshold}',color="gray",
                        alpha=0.5,label=f'N={len(non_group_copy)}',stat="density",kde=True)
            sns.histplot(ax=axes[m, n], data=sub_group_copy, x=f'ratio_filtered_std2_rep1_DNA_{threshold}', color="red",
                            label=f'N={len(sub_group_copy)}',stat="density",kde=True)
            axes[m, n].legend(fontsize="2")
            axes[-1, n].set_xlabel(f'{activity_level}',fontsize=10)
            axes[m, 0].set_ylabel(f'{threshold}',fontsize=10)
            #axes[-1, n].set(xlabel=f'{activity_level}')
            #axes[m, 0].set(ylabel=f'{threshold}')

    fig.suptitle(f'Rep 1 RNA DNA ratio\n{val}', fontsize=12)
    fig.supxlabel('Activity threshold',fontsize=7)
    fig.supylabel('Min DNA count',fontsize=7)
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA_ratio/RNA_DNA_ratio_dist_{val}.png', dpi=1000)

    plt.clf()

################################### RNA DNA CPM #################################################

newpath = f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
    
newpath = f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA/CPM/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
    
 # Loop over keys in the oligo_groups dictionary and create RNA vs DNA graphs   
for key,val in oligo_groups.items():
    fig, axes = plt.subplots(4, 3, sharey=True, sharex=True)
    # Filter the DataFrame into two parts: one for the selected rows and one for the rest
    sub_group = ratio_df[ratio_df['oligo'].str.contains(key, na=False)].copy()
    non_group = ratio_df[~ratio_df['oligo'].str.contains(key, na=False)].copy()
    if key == "a3" or key == "a4":
        sub_group = ratio_df[
            ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        non_group = ratio_df[
            ~ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
    # Loop through the activity levels and thresholds to create subplots
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                non_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = non_group[f"ratio_log_filtered_std2_{rep}"].where(non_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA)
                sub_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = sub_group[f"ratio_log_filtered_std2_{rep}"].where(sub_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA)
                
                if activity_level=="none":
                    non_group_copy = non_group.copy()
                    sub_group_copy = sub_group.copy()
                else:
                    non_group_copy = non_group[(non_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(non_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(non_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    sub_group_copy = sub_group[(sub_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(sub_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(sub_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    
            sns.scatterplot(ax=axes[m, n], data=non_group_copy, x=f'DNA_filtered_std2_cpm_rep1', y=f'RNA_filtered_std2_cpm_rep1', s=5, color="gray", alpha=0.5,label=f'N={len(non_group_copy)}')
            sns.scatterplot(ax=axes[m, n], data=sub_group_copy, x=f'DNA_filtered_std2_cpm_rep1', y=f'RNA_filtered_std2_cpm_rep1', s=5, color="red",label=f'N={len(sub_group_copy)}')            
            axes[m, n].legend(fontsize="2")
            axes[-1, n].set_xlabel(f'{activity_level}',fontsize=10)
            axes[m, 0].set_ylabel(f'{threshold}',fontsize=10)
            # Get axis limits
            min_val = np.nanmin(
                [non_group_copy[f'DNA_filtered_std2_cpm_rep1'].min(), sub_group_copy[f'DNA_filtered_std2_cpm_rep1'].min()]
            )
            max_val = np.nanmax(
                [non_group_copy[f'RNA_filtered_std2_cpm_rep1'].max(), sub_group_copy[f'RNA_filtered_std2_cpm_rep1'].max()]
            )

            # Plot the y = x line
            axes[m, n].plot(
                [min_val, max_val], [min_val, max_val],
                ls='--', color='black', linewidth=1
            )
            #axes[-1, n].set(xlabel=f'{activity_level}')
            #axes[m, 0].set(ylabel=f'{threshold}')

    fig.suptitle(f'RNA vs DNA (CPM) for rep 1\n{val}', fontsize=12)
    fig.supxlabel('Activity threshold', fontsize=7)
    fig.supylabel('Min DNA count',fontsize=7)
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA/CPM/RNA_vs_DNA_filtered_{val}.png', dpi=1000)

    plt.clf()
    
    ################################### RNA DNA raw #################################################

newpath = f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
    
newpath = f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA/raw/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
    
 # Loop over keys in the oligo_groups dictionary and create RNA vs DNA graphs   
for key,val in oligo_groups.items():
    fig, axes = plt.subplots(4, 3, sharey=True, sharex=True)
    # Filter the DataFrame into two parts: one for the selected rows and one for the rest
    sub_group = ratio_df[ratio_df['oligo'].str.contains(key, na=False)].copy()
    non_group = ratio_df[~ratio_df['oligo'].str.contains(key, na=False)].copy()
    if key == "a3" or key == "a4":
        sub_group = ratio_df[
            ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        non_group = ratio_df[
            ~ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
    # Loop through the activity levels and thresholds to create subplots
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                non_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = non_group[f"ratio_log_filtered_std2_{rep}"].where(non_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA)
                sub_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = sub_group[f"ratio_log_filtered_std2_{rep}"].where(sub_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold, pd.NA)
                
                if activity_level=="none":
                    non_group_copy = non_group.copy()
                    sub_group_copy = sub_group.copy()
                else:
                    non_group_copy = non_group[(non_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(non_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(non_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    sub_group_copy = sub_group[(sub_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) |(sub_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) |(sub_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    
            sns.scatterplot(ax=axes[m, n], data=non_group_copy, x=f'DNA_filtered_std2_sum_rep1', y=f'RNA_filtered_std2_sum_rep1', s=5, color="gray", alpha=0.5,label=f'N={len(non_group_copy)}')
            sns.scatterplot(ax=axes[m, n], data=sub_group_copy, x=f'DNA_filtered_std2_sum_rep1', y=f'RNA_filtered_std2_sum_rep1', s=5, color="red",label=f'N={len(sub_group_copy)}')            
            axes[m, n].legend(fontsize="2")
            axes[-1, n].set_xlabel(f'{activity_level}',fontsize=10)
            axes[m, 0].set_ylabel(f'{threshold}',fontsize=10)
            # Get axis limits
            min_val = np.nanmin(
                [non_group_copy[f'DNA_filtered_std2_cpm_rep1'].min(),
                 sub_group_copy[f'DNA_filtered_std2_cpm_rep1'].min()]
            )
            max_val = np.nanmax(
                [non_group_copy[f'RNA_filtered_std2_cpm_rep1'].max(),
                 sub_group_copy[f'RNA_filtered_std2_cpm_rep1'].max()]
            )

            # Plot the y = x line
            axes[m, n].plot(
                [min_val, max_val], [min_val, max_val],
                ls='--', color='black', linewidth=1
            )
            #axes[-1, n].set(xlabel=f'{activity_level}')
            #axes[m, 0].set(ylabel=f'{threshold}')

    fig.suptitle(f'RNA vs DNA (raw) for rep 1\n{val}', fontsize=12)
    fig.supxlabel('Activity threshold', fontsize=7)
    fig.supylabel('Min DNA count',fontsize=7)
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/RNA_DNA/raw/RNA_vs_DNA_filtered_{val}.png', dpi=1000)

    plt.clf()
    
################################### DNA DNA corr #################################################
    
    
newpath = f'./{cells}/{library}{adaptor}/output/filter/DNA_DNA/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
# Rep1 DNA cpm vs Rep2 DNA cpm
for key, val in oligo_groups.items():
    fig, axes = plt.subplots(4, 3, sharey=True, sharex=True)
    # Filter the DataFrame into two parts: one for the selected rows and one for the rest
    sub_group = ratio_df[ratio_df['oligo'].str.contains(key, na=False)].copy()
    non_group = ratio_df[~ratio_df['oligo'].str.contains(key, na=False)].copy()
    if key == "a3" or key == "a4":
        sub_group = ratio_df[
            ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        non_group = ratio_df[
            ~ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
    # Loop through the activity levels and thresholds to create subplots
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                non_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = non_group[
                    f"ratio_log_filtered_std2_{rep}"].where(non_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold,
                                                            pd.NA)
                sub_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = sub_group[
                    f"ratio_log_filtered_std2_{rep}"].where(sub_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold,
                                                            pd.NA)

                if activity_level == "none":
                    non_group_copy = non_group.copy()
                    sub_group_copy = sub_group.copy()
                else:
                    non_group_copy = non_group[
                        (non_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) | (
                                    non_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) | (
                                    non_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    sub_group_copy = sub_group[
                        (sub_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) | (
                                    sub_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) | (
                                    sub_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()

            sns.scatterplot(ax=axes[m, n], data=non_group_copy, x=f'DNA_filtered_std2_cpm_rep1',
                            y=f'DNA_filtered_std2_cpm_rep2', s=5, color="gray", alpha=0.5,
                            label=f'N={len(non_group_copy)}')
            sns.scatterplot(ax=axes[m, n], data=sub_group_copy, x=f'DNA_filtered_std2_cpm_rep1',
                            y=f'DNA_filtered_std2_cpm_rep2', s=5, color="red", label=f'N={len(sub_group_copy)}')
            axes[m, n].legend(fontsize="2")
            axes[-1, n].set_xlabel(f'{activity_level}',fontsize=10)
            axes[m, 0].set_ylabel(f'{threshold}',fontsize=10)
            #axes[-1, n].set(xlabel=f'{activity_level}')
            #axes[m, 0].set(ylabel=f'{threshold}')

    fig.suptitle(f'Correlation Between reps 1 and 2 DNA cpm\n{val}', fontsize=12)
    fig.supxlabel('Activity threshold',fontsize=7)
    fig.supylabel('Min DNA count',fontsize=7)
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/DNA_DNA/DNA_vs_DNA_filtered_{val}.png', dpi=1000)

    plt.clf()
newpath = f'./{cells}/{library}{adaptor}/output/filter/RNA_RNA/'
if not os.path.exists(newpath):
    os.makedirs(newpath)
# Rep1 RNA cpm vs Rep2 RNA cpm
for key, val in oligo_groups.items():
    fig, axes = plt.subplots(4, 3, sharey=True, sharex=True)
    # Filter the DataFrame into two parts: one for the selected rows and one for the rest
    sub_group = ratio_df[ratio_df['oligo'].str.contains(key, na=False)].copy()
    non_group = ratio_df[~ratio_df['oligo'].str.contains(key, na=False)].copy()
    if key == "a3" or key == "a4":
        sub_group = ratio_df[
            ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
        non_group = ratio_df[
            ~ratio_df['oligo'].str.contains(key, na=False) & ~ratio_df['oligo'].str.contains('satmutPilot')].copy()
    # Loop through the activity levels and thresholds to create subplots
    for n, activity_level in enumerate(["none", 0, 1]):
        for m, threshold in enumerate([0, 5, 10, 20]):
            for rep in ['rep1', 'rep2', 'rep3']:
                non_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = non_group[
                    f"ratio_log_filtered_std2_{rep}"].where(non_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold,
                                                            pd.NA)
                sub_group[f'ratio_filtered_std2_{rep}_DNA_{threshold}'] = sub_group[
                    f"ratio_log_filtered_std2_{rep}"].where(sub_group[f"DNA_filtered_std2_sum_{rep}"] >= threshold,
                                                            pd.NA)

                if activity_level == "none":
                    non_group_copy = non_group.copy()
                    sub_group_copy = sub_group.copy()
                else:
                    non_group_copy = non_group[
                        (non_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) | (
                                    non_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) | (
                                    non_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()
                    sub_group_copy = sub_group[
                        (sub_group[f'ratio_filtered_std2_rep1_DNA_{threshold}'] > activity_level) | (
                                    sub_group[f'ratio_filtered_std2_rep2_DNA_{threshold}'] > activity_level) | (
                                    sub_group[f'ratio_filtered_std2_rep3_DNA_{threshold}'] > activity_level)].copy()

            sns.scatterplot(ax=axes[m, n], data=non_group_copy, x=f'RNA_filtered_std2_cpm_rep1',
                            y=f'RNA_filtered_std2_cpm_rep2', s=5, color="gray", alpha=0.5,
                            label=f'N={len(non_group_copy)}')
            sns.scatterplot(ax=axes[m, n], data=sub_group_copy, x=f'RNA_filtered_std2_cpm_rep1',
                            y=f'RNA_filtered_std2_cpm_rep2', s=5, color="red", label=f'N={len(sub_group_copy)}')
            axes[m, n].legend(fontsize="3")
            axes[-1, n].set_xlabel(f'{activity_level}',fontsize=10)
            axes[m, 0].set_ylabel(f'{threshold}',fontsize=10)
            #axes[-1, n].set(xlabel=f'{activity_level}')
            #axes[m, 0].set(ylabel=f'{threshold}')

    fig.suptitle(f'Correlation Between reps 1 and 2 RNA cpm\n{val}', fontsize=12)
    fig.supxlabel('Activity threshold',fontsize=7)
    fig.supylabel('Min DNA count',fontsize=7)
    plt.savefig(f'./{cells}/{library}{adaptor}/output/filter/RNA_RNA/RNA_vs_RNA_filtered_{val}.png', dpi=1000)

    plt.clf()

