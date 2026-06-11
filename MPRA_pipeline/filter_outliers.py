import pandas as pd
import numpy as np
import sys
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import os
import seaborn as sns
from collections import Counter

# have to remove also corresponding barocde from oligo_bs and bc for explode to work properly

def min_max_outliers(rna, dna, oligo_bc):
    # Find indexes corresponding to max and min values
    max_index = rna.index(max(rna))
    min_index = rna.index(min(rna))
    # Remove max and min values from RNA and DNA lists
    rna_filter = [i for j, i in enumerate(rna) if j not in [min_index, max_index]]
    dna_filter = [i for j, i in enumerate(dna) if j not in [min_index, max_index]]
    oligo_bc_filter = [i for j, i in enumerate(oligo_bc) if j not in [min_index, max_index]]
    return rna_filter, dna_filter, oligo_bc_filter

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

# manually copied read counts form alignment output
def create_ratio_df():
    if cells == "chondrocytes" and library == "L1" and adaptor =="a2":
        #chondro L1a2
        rep1DNA = 102933421
        rep1RNA = 211102788
        rep2DNA = 99608813 
        rep2RNA = 252711591
        rep3DNA = 96024029 
        rep3RNA = 285094704
        
    if cells == "chondrocytes" and library == "L2" and adaptor =="a2":
    #   chondro L2a2
        rep1DNA = 70575053
        rep1RNA = 237195340 
        rep2DNA = 88628404
        rep2RNA = 252883941
        rep3DNA = 79947531
        rep3RNA = 333520513

    if cells == "chondrocytes" and library == "L3" and adaptor =="a2":
        #chondro L3a2
        rep1DNA = 86012916
        rep1RNA = 262288284
        rep2DNA = 84852307
        rep2RNA = 267028259
        rep3DNA = 64816951
        rep3RNA = 252832900
        
    if cells == "chondrocytes" and library == "L1" and adaptor =="a1":
        #chondro L1a1
        rep1DNA = 81111655
        rep1RNA = 251065756
        rep2DNA = 82290336 
        rep2RNA = 281299871
        rep3DNA = 78159990 
        rep3RNA = 305537255
        
    if cells == "chondrocytes" and library == "L2" and adaptor =="a1":
        #chondro L1a2
        rep1DNA = 100863191
        rep1RNA = 279158189
        rep2DNA = 78870276 
        rep2RNA = 288849225
        rep3DNA = 73624204 
        rep3RNA = 288161895

    if cells == "neurons" and library == "L1" and adaptor =="a1":
        #neurons L1a1
        rep1DNA_first = 43955698
        rep1RNA_first = 136943218
        rep2DNA_first = 47316452
        rep2RNA_first = 141580219
        rep3DNA_first = 46544258
        rep3RNA_first = 145645192
        rep1DNA_second = 74092786 
        rep1RNA_second = 204129993 
        rep2DNA_second = 77150709 
        rep2RNA_second = 206148983 
        rep3DNA_second = 73081201 
        rep3RNA_second = 249183776

        rep1DNA = rep1DNA_first + rep1DNA_second
        rep1RNA = rep1RNA_first + rep1RNA_second
        rep2DNA = rep2DNA_first + rep2DNA_second
        rep2RNA = rep2RNA_first + rep2RNA_second
        rep3DNA = rep3DNA_first + rep3DNA_second
        rep3RNA = rep3RNA_first + rep3RNA_second
        
    ### 

    # access replicate UMI counts files
    rep1 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep1.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])
    rep2 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep2.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])
    rep3 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep3.txt', sep='\t', header=0, names=["oligo_bc_no_filter", "RNA_no_filter", "DNA_no_filter"])



    # get number of reads reads
    DNA_reads = [rep1DNA, rep2DNA, rep3DNA]
    RNA_reads = [rep1RNA, rep2RNA, rep3RNA]
    read_dict = {"DNA_reads": DNA_reads, "RNA_reads": RNA_reads}
    print(DNA_reads)
    print(RNA_reads)

    # group oligos and sum RNA, DNA counts. Also counts how many barcodes we have per oligo. Calculate ratio and perform log2 on ratio
    replicates = []
    exploded_df = []
    
    for i, rep in enumerate([rep1, rep2, rep3]):
        for outlier_filter in ["filtered_min_max", "filtered_std2", "filtered_std3"]:
            rep[f'oligo_bc_{outlier_filter}'] = None
            for molecule in ["RNA", "DNA"]:
                rep[f'{molecule}_{outlier_filter}'] = None
        rep = rep[(rep['DNA_no_filter'] > 0) | (rep['RNA_no_filter'] > 0)]
        rep[['oligo', 'bc']] = rep.oligo_bc_no_filter.str.rsplit("_", 1, expand=True)
        rep = rep.set_index('oligo')
        rep_list = rep.groupby('oligo').aggregate(lambda a: a.tolist())
        rep_list = rep_list.reset_index()
        for index, row in rep_list.iterrows():
            rep_list.at[index, 'RNA_filtered_min_max'], rep_list.at[index, 'DNA_filtered_min_max'], rep_list.at[index, 'oligo_bc_filtered_min_max'] = min_max_outliers(
                row['RNA_no_filter'], row['DNA_no_filter'], row['oligo_bc_no_filter'])
            rep_list.at[index, 'RNA_filtered_std2'], rep_list.at[index, 'DNA_filtered_std2'], rep_list.at[index, 'oligo_bc_filtered_std2'] = std_n_outliers(
                row['RNA_no_filter'], row['DNA_no_filter'], row['oligo_bc_no_filter'], 2)
            rep_list.at[index, 'RNA_filtered_std3'], rep_list.at[index, 'DNA_filtered_std3'], rep_list.at[index, 'oligo_bc_filtered_std3'] = std_n_outliers(
                row['RNA_no_filter'], row['DNA_no_filter'], row['oligo_bc_no_filter'], 3)

        for outlier_filter in ["no_filter", "filtered_min_max", "filtered_std2", "filtered_std3"]:
            exploded = rep_list.explode([f"RNA_{outlier_filter}", f"DNA_{outlier_filter}", f'oligo_bc_{outlier_filter}'])
            exploded_df.append(exploded)
            rep_list = rep_list.drop([f'oligo_bc_{outlier_filter}'], axis=1)
        rep_list = rep_list.drop(['bc'], axis=1)
        for outlier_filter in ["no_filter", "filtered_min_max", "filtered_std2", "filtered_std3"]:
            rep_list[f'count_{outlier_filter}'] = rep_list[f'RNA_{outlier_filter}'].apply(len)
            mask = (rep_list[f'count_{outlier_filter}'] > 0) # filtering might lead to some cases where a barcode ends up with 0 barcodes - want to filter out these oligos
            rep_valid = rep_list[mask]
            for molecule in ["RNA", "DNA"]:
                rep_list.loc[mask, f'{molecule}_{outlier_filter}_sum'] = rep_list.loc[mask, f'{molecule}_{outlier_filter}'].apply(sum)
                rep_list.loc[mask, f"{molecule}_{outlier_filter}_pseudo"] = rep_list.loc[mask, f"{molecule}_{outlier_filter}_sum"] + 1
                rep_list.loc[mask, f"{molecule}_{outlier_filter}_cpm"] = (rep_list.loc[mask, f"{molecule}_{outlier_filter}_pseudo"] * 1000000)/read_dict[f"{molecule}_reads"][0] # this has to be fixed - change 0 to i

            rep_list.loc[mask, f"ratio_{outlier_filter}"] = rep_list.loc[mask, f"RNA_{outlier_filter}_cpm"] / rep_list.loc[mask, f"DNA_{outlier_filter}_cpm"]
            rep_list.loc[mask, f"ratio_{outlier_filter}_log"] = np.log2(rep_list.loc[mask, f"ratio_{outlier_filter}"])
        replicates.append(rep_list)

    rep1_counts = replicates[0]
    rep2_counts = replicates[1]
    rep3_counts = replicates[2]

    # merge the different outputs
    ratio_df = rep1_counts.merge(rep2_counts,on='oligo',how="outer",suffixes=("_rep1",None)).merge(rep3_counts,on='oligo',how="outer",suffixes=("_rep2","_rep3"))
    # save df
    ratio_df.to_csv(f'./{cells}/{library}{adaptor}/output/activity/ratio_wo_outliers.csv', header=True, index=False)
    
    # todo: save the exploded df (2std) - so we get a file in the same format as barcode_counts_UMI.txt - which than can be used to create the inout for MPRAnalyze and the other steps of the analysis

    return ratio_df

if not os.path.exists(f'./{cells}/{library}{adaptor}/output/activity/ratio_wo_outliers.csv'):
    ratio_df = create_ratio_df()
else:
    ratio_df = pd.read_csv(f'./{cells}/{library}{adaptor}/output/activity/ratio_wo_outliers.csv')
    
    
# # overview over how many we loose per oligo for the different methods
# plt.clf()
# fig, axes = plt.subplots(3, 3, sharey=True)
# for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    # for m, outlier_filter in enumerate(["filtered_min_max","filtered_std2","filtered_std3"]):
        # axes[n,m].hist(ratio_df[f"count_no_filter_{rep}"]-ratio_df[f"count_{outlier_filter}_{rep}"], bins=25, range=(0,24))
        # print(Counter(ratio_df[f"count_no_filter_{rep}"]-ratio_df[f"count_{outlier_filter}_{rep}"]))
# plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/histogram_number_filter.pdf')

# # calculate correlation between different methods within rep
# print("corr between methods")
# plt.clf()
# fig, axes = plt.subplots(1,3, sharey=True)
# cbar_ax = fig.add_axes([.91, .3, .03, .4])
# x_label_list = ['no_filter', 'min_max', 'std2', 'std3']
# y_label_list = ['no_filter', 'min_max', 'std2', 'std3']
# for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    # corr = ratio_df[[f"ratio_no_filter_log_{rep}",f"ratio_filtered_min_max_log_{rep}",f"ratio_filtered_std2_log_{rep}",f"ratio_filtered_std3_log_{rep}"]].corr()
    # print(corr.to_string())
    # sns.heatmap(ax=axes[n], data=corr, cmap="Reds", vmin=0, vmax=1, annot=True, annot_kws={"fontsize":8}, fmt=".2f", square=True, xticklabels=x_label_list, yticklabels=y_label_list, cbar_ax=cbar_ax)
    # axes[n].tick_params(axis='x', labelrotation = 90)
# plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/correlation_between_methods_of_outlier_removal.pdf')


# # histogram of barcode counts and of dna counts
# string_list=["", "sum_"]
# for xlim in ["", "xlim20"]:
    # for i, counts in enumerate(["count", "DNA"]):
        # plt.clf()
        # fig, axes = plt.subplots(3, 4, sharey=True)
        # for n, rep in enumerate(["rep1", "rep2", "rep3"]):
            # for m, outlier_filter in enumerate(["no_filter", "filtered_min_max","filtered_std2","filtered_std3"]):
                # if xlim == "xlim20":
                    # x, bins, patches = axes[n,m].hist(ratio_df[f"{counts}_{outlier_filter}_{string_list[i]}{rep}"], range = (0,19), bins=20)
                # else:
                    # x, bins, patches = axes[n,m].hist(ratio_df[f"{counts}_{outlier_filter}_{string_list[i]}{rep}"])
                # print(bins)
        # plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/histogram_{counts}_filter_{xlim}.pdf')


# correlation between the reps for the different methods and different bc/dns count thresholds and activity thresholds
string_list=["", "sum_"]
for activity_level in ["none", 0, 1]:
    for i, counts in enumerate(["count", "DNA"]):
        plt.clf()
        fig, axes = plt.subplots(4,4, sharey=True, sharex=True)
        cbar_ax = fig.add_axes([.91, .3, .03, .4])
        x_label_list = ['rep1', 'rep2', 'rep3']
        y_label_list = ['rep1', 'rep2', 'rep3']
        for n, outlier_filter in enumerate(["no_filter","filtered_min_max","filtered_std2","filtered_std3"]):
            for m, threshold in enumerate([0, 5, 10, 20]):
                for rep in ['rep1', 'rep2', 'rep3']:
                    ratio_df[f'ratio_{outlier_filter}_{rep}_{counts}_{threshold}'] = ratio_df[f"ratio_{outlier_filter}_log_{rep}"].where(ratio_df[f"{counts}_{outlier_filter}_{string_list[i]}{rep}"] >= threshold, pd.NA)
                if activity_level=="none":
                    ratio_df_copy = ratio_df.copy()
                else:
                    ratio_df_copy = ratio_df[(ratio_df[f'ratio_{outlier_filter}_rep1_{counts}_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep2_{counts}_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep3_{counts}_{threshold}'] > activity_level)].copy()
                corr = ratio_df_copy[[f'ratio_{outlier_filter}_rep1_{counts}_{threshold}',f'ratio_{outlier_filter}_rep2_{counts}_{threshold}',f'ratio_{outlier_filter}_rep3_{counts}_{threshold}']].corr()
                sns.heatmap(ax=axes[m,n], data=corr, cmap="Reds", vmin=0, vmax=1, annot=True, annot_kws={"fontsize":6}, fmt=".2f", square=True, xticklabels=x_label_list, yticklabels=y_label_list, cbar_ax=cbar_ax)
                axes[m, n].tick_params(axis='x', labelrotation = 90)
           
        # plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/correlation_between_reps_with_outlier_removal_and_{counts}_filter_activityfilter={activity_level}.pdf')
        plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/correlation_between_reps_with_outlier_removal_and_{counts}_filter_activityfilter={activity_level}.png', dpi=1000)

# scatter plots of correlation between rep1 and rep2 for the different methods and different bc/dns count thresholds and activity thresholds
string_list=["", "sum_"]
for activity_level in ["none", 0, 1]:
    for i, counts in enumerate(["count", "DNA"]):
        plt.clf()
        fig, axes = plt.subplots(4,4, sharey=True, sharex=True)
        for n, outlier_filter in enumerate(["no_filter","filtered_min_max","filtered_std2","filtered_std3"]):
            for m, threshold in enumerate([0, 5, 10, 20]):
                for rep in ['rep1', 'rep2', 'rep3']:
                    ratio_df[f'ratio_{outlier_filter}_{rep}_{counts}_{threshold}'] = ratio_df[f"ratio_{outlier_filter}_log_{rep}"].where(ratio_df[f"{counts}_{outlier_filter}_{string_list[i]}{rep}"] >= threshold, pd.NA)
                if activity_level=="none":
                    ratio_df_copy = ratio_df.copy()
                else:
                    ratio_df_copy = ratio_df[(ratio_df[f'ratio_{outlier_filter}_rep1_{counts}_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep2_{counts}_{threshold}'] > activity_level) |(ratio_df[f'ratio_{outlier_filter}_rep3_{counts}_{threshold}'] > activity_level)].copy()
                sns.scatterplot(ax=axes[m, n], data = ratio_df_copy, x = f'ratio_{outlier_filter}_rep1_{counts}_{threshold}', y = f'ratio_{outlier_filter}_rep2_{counts}_{threshold}', s=5)
                print(counts)
                print(outlier_filter)
                print(threshold)
           
        # plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/scatter_correlation_between_reps_with_outlier_removal_and_{counts}_filter_activityfilter={activity_level}.pdf')
        plt.savefig(f'./{cells}/{library}{adaptor}/output/activity/scatter_correlation_between_reps_with_outlier_removal_and_{counts}_filter_activityfilter={activity_level}.png', dpi=1000)

# # calculate oligo coverage
# print("oligo coverage")
# string_list=["", "sum_"]
# for i, counts in enumerate(["count", "DNA"]):
    # print(counts)
    # for n, outlier_filter in enumerate(["no_filter","filtered_min_max","filtered_std2","filtered_std3"]):
        # print(outlier_filter)
        # for m, threshold in enumerate([0, 5, 10, 20]):
            # print(threshold)
            # for rep in ["rep1", "rep2", "rep3"]:
                # count = (ratio_df[f"{counts}_{outlier_filter}_{string_list[i]}{rep}"] >= threshold).sum()
                # ratio = round(count/79812, 2)
                # print(rep, ratio)
                