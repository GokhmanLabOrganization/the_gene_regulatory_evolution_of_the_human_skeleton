# Creating file with quantitative results

# what to include in df:
## oligo name x
## rna/dna ratio for replicates combined - how to combine replicates? sum reads? x
## rna/dna ratio for each replicate alone x
## activity based on alpha x
## control annotation x
## output from mpranalyze quantitative x 
## barcode statistics - what exactly? - possible barcodes, RNA/DNA barcodes in each rep,  x
## UMI statistics - what exactly - list with UMI counts for each barcode in each replicate, maxUMI? x

# figures to create:
## corr between replicates rna/dna
## corr only active
## corr between pairs - both for alpha and rna/dna?
## controls boxplot - both for alpha and rna/dna?
## alpha distribution per group (boxplot)
## distribution of alpha
## distribution of rna/dna
## distribution of rna counts (normalized to number of barcodes)
## distribution of dna counts (normalized to number of barcodes)
## alpha vs rna/dna
## UMIs per barcodes plots
## barcodes/oligo statistics graph
## alpha vs # bcs???
## alpha vs max UMI count???

##TODO
## 10/06/24: some analyses use pseudo counts but they do not exist in the dataframe anymore (solution: change the analyses to use the raw counts)

# question - for rna/dna - all graphs for each rep separately and for all reps combined?


# other outputs:
## oligo coverage
## barcodes/oligo statistics
## % active

import pandas as pd
import sys
import os
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy import stats
import numpy as np
import re
import statsmodels.stats.multitest as smm
from scipy.stats import pearsonr

def get_fig_xlim_ylim(x_data,y_data):
    curr_minima = np.floor(min(min(x_data),min(y_data)))
    curr_maxima = np.ceil(max(max(x_data), max(y_data)))
    return curr_minima,curr_maxima

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]


activity_folder_path = f'./{cells}/{library}{adaptor}/output/activity_after_filter'
comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv'
#comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr_L3a3_Hh.csv' # NM 01-11-2024 Toggle on to include only Hh oligos for Hh analysis. Toggle off previous line

mpranalyze_output = f'./{cells}/{library}{adaptor}/output/mpranalyze_quantitative/{cells}_{library}{adaptor}_results_full_quantitative_nobc_fdr_annot_factor_new_design.txt'

def create_comb_df():
    # read in mpranalyze quantitative output file
    output = pd.read_csv(mpranalyze_output, sep='\t', header =0)
    output['oligo'] = output.index
    print(output.head().to_string())

    # read in UMI counts file
    # UMI_df = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI.txt', sep='\t')
    UMI_df = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t')

    # read in file with RNA/DNA ratio accross all replicates
    RNA_DNA = pd.read_csv(f'./{cells}/{library}{adaptor}/output/activity_after_filter/ratio_after_filter.csv', usecols=["oligo_column", "RNA_rep1", "DNA_rep1", "RNA_rep2", "DNA_rep2", "RNA_rep3", "DNA_rep3", "RNA_rep_comb", "DNA_rep_comb", "ratio_log_rep1", "ratio_log_rep2", "ratio_log_rep3", "ratio_log_rep_comb", "count_rep1", "count_rep2", "count_rep3", "count_rep_comb"])
    # RNA_DNA = pd.read_csv(f'./{cells}/{library}{adaptor}/output/activity/ratio.csv', usecols=["oligo_column", "RNA_rep1", "DNA_rep1", "RNA_rep2", "DNA_rep2", "RNA_rep3", "DNA_rep3", "ratio_log_rep1", "ratio_log_rep2", "ratio_log_rep3", "ratio_log_rep_comb", "count_rep1", "count_rep2", "count_rep3", "count_rep_comb"])

    # create list of UMI counts
    UMI_df[['oligo','bc']] = UMI_df.oligo_bc.str.rsplit("_", n=1, expand=True) #06-10-2024 N.M Refining the split function behavior by specifying the argument n = 1
    UMI_df = UMI_df.set_index('oligo')
    UMI_df_list = UMI_df.groupby('oligo').aggregate(lambda a: a.tolist())
    UMI_df_list = UMI_df_list.reset_index()
    UMI_df_list = UMI_df_list.drop(['oligo_bc', 'bc'], axis=1)
    print(UMI_df_list.head())


    # merge mprAnalyze quantitiative output file to UMI list
    merged_df = UMI_df_list.merge(output, on = "oligo", how = "left")
    merged_df['possible_bcs'] = merged_df[f'RNA_{cells}_rep1'].str.len()

    # add bcs count for each rep and molecule
    merged_df['bcs_RNA_rep1'] = merged_df[f'RNA_{cells}_rep1'].apply(lambda x: len([a for a in x if a > 0]))
    merged_df['bcs_RNA_rep2'] = merged_df[f'RNA_{cells}_rep2'].apply(lambda x: len([a for a in x if a > 0]))
    merged_df['bcs_RNA_rep3'] = merged_df[f'RNA_{cells}_rep3'].apply(lambda x: len([a for a in x if a > 0]))
    merged_df['bcs_DNA_rep1'] = merged_df[f'DNA_{cells}_rep1'].apply(lambda x: len([a for a in x if a > 0]))
    merged_df['bcs_DNA_rep2'] = merged_df[f'DNA_{cells}_rep2'].apply(lambda x: len([a for a in x if a > 0]))
    merged_df['bcs_DNA_rep3'] = merged_df[f'DNA_{cells}_rep3'].apply(lambda x: len([a for a in x if a > 0]))
    
    merged_df['bcs_RNA_comb']=merged_df[[f'RNA_{cells}_rep1',f'RNA_{cells}_rep2',f'RNA_{cells}_rep3']].apply(lambda x: sum([1 if (x[0][i]!=0 or x[1][i]!=0 or x[2][i]!=0) else 0 for i in range(len(x[0]))])  ,axis=1)
    merged_df['bcs_DNA_comb']=merged_df[[f'DNA_{cells}_rep1',f'DNA_{cells}_rep2',f'DNA_{cells}_rep3']].apply(lambda x: sum([1 if (x[0][i]!=0 or x[1][i]!=0 or x[2][i]!=0) else 0 for i in range(len(x[0]))])  ,axis=1)

    # merge rna/dna ratio to dataframe
    merged_with_ratio = merged_df.merge(RNA_DNA, left_on = "oligo", right_on = "oligo_column", how = "left")

    # add control annotations
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('PosCtrl_osteoblast_active'), 'control_type'] = 'PosCtrl_osteoblast_active'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('PosCtrl_chondrocyte_active'), 'control_type'] = 'PosCtrl_chondrocyte_active'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('PosCtrl_neuron_active'), 'control_type'] = 'PosCtrl_neuron_active'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('NegCtrl_non_SCREEN'), 'control_type'] = 'NegCtrl_non_SCREEN'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('NegCtrl_active_not_diff'), 'control_type'] = 'NegCtrl_active_not_diff'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('scrambled_control'), 'control_type'] = 'scrambled_control'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('NegCtrl_not_active'), 'control_type'] = 'NegCtrl_not_active'
    merged_with_ratio.loc[merged_with_ratio['oligo'].str.contains('PosCtrl_diff'), 'control_type'] = 'PosCtrl_diff'
    merged_with_ratio['control_type'] = merged_with_ratio['control_type'].fillna(value='No control')

    # drop oligo column
    merged_with_ratio = merged_with_ratio.drop(['oligo_column'], axis=1)

    # add active/non active
    merged_with_ratio["activity"] = "non_active"
    merged_with_ratio.loc[merged_with_ratio["fdr.mad"] < 0.05, "activity"] = "active"

    # add count of bc with highest count
    def get_max(col1, col2, col3):
        return max([max(col1), max(col2), max(col3)])

    merged_with_ratio['max_UMI_RNA'] = merged_with_ratio.apply(lambda x: get_max(x[f'RNA_{cells}_rep1'], x[f'RNA_{cells}_rep2'], x[f'RNA_{cells}_rep3']), axis=1)
    merged_with_ratio['max_UMI_DNA'] = merged_with_ratio.apply(lambda x: get_max(x[f'DNA_{cells}_rep1'], x[f'DNA_{cells}_rep2'], x[f'DNA_{cells}_rep3']), axis=1)
    
    # add count of bc with highest count per replicate
    for rep in ["rep1","rep2","rep3"]:
        merged_with_ratio[f'max_UMI_RNA_{rep}'] = merged_with_ratio[f'RNA_{cells}_{rep}'].apply(max)
        merged_with_ratio[f'max_UMI_DNA_{rep}'] = merged_with_ratio[f'DNA_{cells}_{rep}'].apply(max)
        
    # add fdr that does not consider (1) controls, (2) oligos with less than 5 DNA counts (3) L4 unrelated oligos FABP7, hh.SCREEN
    
    filter_mask= ((merged_with_ratio['DNA_rep_comb'] >= 5) &
                    (merged_with_ratio['control_type'] == 'No control')&
                    (~merged_with_ratio['pval.mad'].isna())&
                    (~merged_with_ratio['oligo'].str.contains('_Hh_'))&
                    (~merged_with_ratio['oligo'].str.contains('hh.SCREEN'))& # hh SCREEN elements (without variants)
                    (~merged_with_ratio['oligo'].str.contains('FABP7')))
    
    filtered_merged_with_ratio = merged_with_ratio[filter_mask]
    #filtered_merged_with_ratio = merged_with_ratio[(merged_with_ratio['DNA_rep_comb'] >= 5) & (merged_with_ratio['control_type'] == 'No control')& (~merged_with_ratio['pval.mad'].isna())& (~merged_with_ratio['oligo'].str.contains('_SCREEN_'))] # NM 01-11-2024 Toggle on to include only Hh oligos for Hh analysis. Toggle off previous line
    print(len(filtered_merged_with_ratio.index))
    rej, pval_corr = smm.fdrcorrection(filtered_merged_with_ratio["pval.mad"])
    print(len(pval_corr))
    merged_with_ratio.loc[filter_mask, 'fdr.mad_adjusted'] = pval_corr
    #merged_with_ratio.loc[(merged_with_ratio['DNA_rep_comb'] >= 5) & (merged_with_ratio['control_type'] == 'No control')& (~merged_with_ratio['pval.mad'].isna())& (~merged_with_ratio['oligo'].str.contains('_SCREEN_')), 'fdr.mad_adjusted'] = pval_corr # NM 01-11-2024 Toggle on to include only Hh oligos for Hh analysis. Toggle off previous line

    merged_with_ratio['fdr.mad_adjusted'].fillna(1, inplace=True)
    
    # # add fdr that should be the same as original fdr as test
    
    # test_merged_with_ratio = merged_with_ratio[(~merged_with_ratio['fdr.mad'].isna())]
    # print(len(test_merged_with_ratio.index))
    # rej, pval_corr = smm.fdrcorrection(test_merged_with_ratio["pval.mad"])
    # print(len(pval_corr))
    # merged_with_ratio.loc[(~merged_with_ratio['fdr.mad'].isna()), 'fdr.test'] = pval_corr
    # merged_with_ratio['fdr.test'].fillna(1, inplace=True)
    
    
    # add new activity columns
    merged_with_ratio.loc[merged_with_ratio["fdr.mad_adjusted"] < 0.05, "activity_adjusted"] = "active"
    merged_with_ratio.loc[(merged_with_ratio['DNA_rep_comb'] >= 5) &
                            (merged_with_ratio['control_type'] == 'No control')&
                            (~merged_with_ratio['pval.mad'].isna()) &
                            (merged_with_ratio["fdr.mad_adjusted"] >= 0.05)&
                            (~merged_with_ratio['oligo'].str.contains('hh.SCREEN'))&
                            (~merged_with_ratio['oligo'].str.contains('FABP7')) &
                            (~merged_with_ratio['oligo'].str.contains('_Hh_')), 'activity_adjusted'] = "non_active"
    #merged_with_ratio.loc[(merged_with_ratio['DNA_rep_comb'] >= 5) & (merged_with_ratio['control_type'] == 'No control')& (~merged_with_ratio['pval.mad'].isna()) & (merged_with_ratio["fdr.mad_adjusted"] >= 0.05)& (~merged_with_ratio['oligo'].str.contains('_SCREEN_')), 'activity_adjusted'] = "non_active" # NM 01-11-2024 Toggle on to include only Hh oligos for Hh analysis. Toggle off previous line
    
    # add column whether oligo is input or not for mpranalyze comp
    merged_with_ratio["input_comparative"] = "no"
    merged_with_ratio.loc[merged_with_ratio["fdr.mad_adjusted"] < 0.05, "input_comparative"] = "yes"
    merged_with_ratio.loc[merged_with_ratio['control_type'] == "PosCtrl_diff", "input_comparative"] = "yes"
    
    
    # # qc: add column for difference of fdr and fdr adjusted
    # merged_with_ratio["fdr-fdradjusted"] = merged_with_ratio["fdr.mad"] - merged_with_ratio["fdr.mad_adjusted"]
    # save
    merged_with_ratio.to_csv(comb_df_path, header=True, index = False)
    
    return merged_with_ratio

if not os.path.exists(comb_df_path):
    data_df = create_comb_df()
else:
    data_df = pd.read_csv(comb_df_path)
    
    
# count active/non active (MAIN ANALYSIS)
# still open question: how to calculate percent active properly - mostly what should be the denominator?
num_active = (data_df["activity"] == "active").sum()
num_not_active = (data_df["activity"] == "non_active").sum()
print("before adjusting fdr")
print(f'# active: {num_active}')
print(f'# not active: {num_not_active}')
print(f'% active: {num_active/(num_active+num_not_active)}')

print("after adjusting fdr")
num_active = (data_df["activity_adjusted"] == "active").sum()
num_not_active = (data_df["activity_adjusted"] == "non_active").sum()
print(f'# active: {num_active}')
print(f'# not active: {num_not_active}')
print(f'% active: {num_active/(num_active+num_not_active)}')

# todo: add percent active also after adjusting fdr and filter on number of dna counts

# plots
sns.set(font_scale=0.35)
sns.set_style("whitegrid")

# boxplot - activity of different controls groups(MAIN ANALYSIS)
plt.clf()
sns.boxplot(data=data_df, x="alpha", y="control_type", order=['PosCtrl_neuron_active', 'PosCtrl_osteoblast_active', 'PosCtrl_chondrocyte_active', 'NegCtrl_active_not_diff', 'PosCtrl_diff', 'NegCtrl_not_active', 'scrambled_control', 'NegCtrl_non_SCREEN', 'No control'], palette=["gold", "gold", "gold", "gold", "gold", "silver", "silver", "silver", "lightcoral"])
plt.savefig(f'{activity_folder_path}/alpha_per_control_type_boxplot_{cells}_{library}{adaptor}.png', dpi=330)
plt.savefig(f'{activity_folder_path}/alpha_per_control_type_boxplot_{cells}_{library}{adaptor}.pdf')

plt.clf()
sns.boxplot(data=data_df, x="alpha", y="control_type", showfliers = False, order=['PosCtrl_neuron_active', 'PosCtrl_osteoblast_active', 'PosCtrl_chondrocyte_active', 'NegCtrl_active_not_diff', 'PosCtrl_diff', 'NegCtrl_not_active', 'scrambled_control', 'NegCtrl_non_SCREEN', 'No control'], palette=["gold", "gold", "gold", "gold", "gold", "silver", "silver", "silver", "lightcoral"])
plt.savefig(f'{activity_folder_path}/alpha_per_control_type_boxplot_{cells}_{library}{adaptor}_no_outliers.pdf')
plt.savefig(f'{activity_folder_path}/alpha_per_control_type_boxplot_{cells}_{library}{adaptor}_no_outliers.png', dpi=330)


# # correlation between reps (MAIN ANALYSIS)
plt.clf()
fig, axes = plt.subplots(1,3)
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    print(rep)
    corr_all = data_df.loc[:, f"ratio_log_{rep[0]}"].corr(data_df.loc[:, f"ratio_log_{rep[1]}"])
    corr_5 = data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5), f"ratio_log_{rep[1]}"])
    corr_10 = data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 10)& (data_df[f"DNA_{rep[1]}"] >= 10), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 10)& (data_df[f"DNA_{rep[1]}"] >= 10), f"ratio_log_{rep[1]}"])

    data_df_temp = data_df.dropna(subset = [f"ratio_log_{rep[0]}", f"ratio_log_{rep[1]}"])

    print("all oligos:", corr_all)
    print(pearsonr(data_df_temp.loc[:, f"ratio_log_{rep[0]}"], data_df_temp.loc[:, f"ratio_log_{rep[1]}"]))
    print(">=5DNA:", corr_5)
    print(pearsonr(data_df_temp.loc[(data_df_temp[f"DNA_{rep[0]}"] >= 5)& (data_df_temp[f"DNA_{rep[1]}"] >= 5), f"ratio_log_{rep[0]}"], data_df_temp.loc[(data_df_temp[f"DNA_{rep[0]}"] >= 5)& (data_df_temp[f"DNA_{rep[1]}"] >= 5), f"ratio_log_{rep[1]}"]))
    print(">=10DNA:", corr_10)
    
    values = np.vstack([data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)][f"ratio_log_{rep[0]}"], data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)][f"ratio_log_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    sns.scatterplot(ax=axes[n], data = data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", c=kernel, s=1, linewidth=0)
    axes[n].set(xlabel=f'{rep[0]}', ylabel=f'{rep[1]}')
    axes[n].set_aspect('equal')
    axes[n].yaxis.get_label().set_visible(True)
plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - >=5bcs')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_5_DNA_rna_dna.pdf')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_5_DNA_rna_dna.png', dpi=330)

# correlation for only active sequences (here activity is defined according to the fdr that is done on all seqeunces - activity column and not activity_adjusted) (MAIN ANALYSIS)
plt.clf()
fig, axes = plt.subplots(1,3)
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    print(rep)
    corr_active_all = data_df.loc[(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])
    corr_active_5 = data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])
    corr_active_10 = data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 10)& (data_df[f"DNA_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"DNA_{rep[0]}"] >= 10)& (data_df[f"DNA_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])

    data_df_temp = data_df.dropna(subset = [f"ratio_log_{rep[0]}", f"ratio_log_{rep[1]}"])
    
    print("all DNA:", corr_active_all)
    print(pearsonr(data_df_temp.loc[(data_df_temp["activity"] == "active"), f"ratio_log_{rep[0]}"], data_df_temp.loc[(data_df_temp["activity"] == "active"), f"ratio_log_{rep[1]}"]))
    print(">=5DNA:", corr_active_5)
    print(pearsonr(data_df_temp.loc[(data_df_temp[f"DNA_{rep[0]}"] >= 5)& (data_df_temp[f"DNA_{rep[1]}"] >= 5)&(data_df_temp["activity"] == "active"), f"ratio_log_{rep[0]}"], data_df_temp.loc[(data_df_temp[f"DNA_{rep[0]}"] >= 5)& (data_df_temp[f"DNA_{rep[1]}"] >= 5)&(data_df_temp["activity"] == "active"), f"ratio_log_{rep[1]}"]))
    print(">=10DNA:", corr_active_10)
    
    values = np.vstack([data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"ratio_log_{rep[0]}"], data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"ratio_log_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    sns.scatterplot(ax=axes[n], data = data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", c=kernel, s=1, linewidth=0)
    axes[n].set(xlabel=f'{rep[0]}', ylabel=f'{rep[1]}')
    axes[n].set_aspect('equal')
    axes[n].yaxis.get_label().set_visible(True)

plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - active, >=5bcs')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5_DNA_rna_dna.pdf')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5_DNA_rna_dna.png', dpi=330)

# # correlation of osteoblast active controls and neuron active controls (not main analysis)
# for control_type in ["PosCtrl_osteoblast_active", "PosCtrl_neuron_active"]:
    # print(control_type)
    # plt.clf()
    # fig, axes = plt.subplots(1, 3, subplot_kw=dict(box_aspect=1))

    # fig.suptitle(f'Correlation of {control_type} between replicates')
    # for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
        # sns.scatterplot(ax=axes[n], data = data_df[(data_df[f"control_type"] == f"{control_type}")&(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", s=8)
        # print(rep)
        # corr_active_all = data_df.loc[(data_df[f"control_type"] == f"{control_type}")&(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"control_type"] == f"{control_type}")&(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])
        # corr_active_5 = data_df.loc[(data_df[f"control_type"] == f"{control_type}")&(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"control_type"] == f"{control_type}")&(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"ratio_log_{rep[1]}"])
        # print("all bcs:", corr_active_all)
        # print(">=5bcs:", corr_active_5)
        # # axes[n].axis('equal')
    # plt.tight_layout()
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_{control_type}.pdf')

## The follwing analyses create a scatter plot comparing rna/dna ratios across replicates and then color the scatterplot according to different features (not main analysis)
# #correlation plot of rna/dna - colored according to barcode counts
# plt.clf()
# fig, axes = plt.subplots(3, 3)
# fig.suptitle('Correlation of RNA/DNA between replicates, colored by #bcs')
# for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    # sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="count_rep1", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="count_rep2", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="count_rep3", s=5,alpha = .9)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_bcs_.pdf')

# # correlation plot of rna/dna - colored according to dna counts
# plt.clf()
# fig, axes = plt.subplots(3, 3)
# fig.suptitle('Correlation of RNA/DNA between replicates, colored by dna cpm')
# for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    # sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="DNA_pseudo_rep1", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="DNA_pseudo_rep2", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="DNA_pseudo_rep3", s=5,alpha = .9)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_dna.pdf')

# # correlation plot of rna/dna - colored according to rna counts
# plt.clf()
# fig, axes = plt.subplots(3, 3)
# fig.suptitle('Correlation of RNA/DNA between replicates, colored by rna cpm')
# for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    # sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="RNA_pseudo_rep1", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="RNA_pseudo_rep2", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="RNA_pseudo_rep3", s=5,alpha = .9)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_rna.pdf')

# # correlation plot of rna/dna - colored according to maxUMI dna counts
# plt.clf()
# fig, axes = plt.subplots(3, 3)
# fig.suptitle('Correlation of RNA/DNA between replicates, colored by maxUMIdna')
# for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    # sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_DNA_rep1", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_DNA_rep2", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_DNA_rep3", s=5,alpha = .9)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_maxUMIdna.pdf')

# # correlation plot of rna/dna - colored according to maxUMI rna counts
# plt.clf()
# fig, axes = plt.subplots(3, 3)
# fig.suptitle('Correlation of RNA/DNA between replicates, colored by maxUMIrna')
# for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    # sns.scatterplot(ax=axes[0, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_RNA_rep1", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[1, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_RNA_rep2", s=5,alpha = .9)
    # sns.scatterplot(ax=axes[2, n], data = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", hue="max_UMI_RNA_rep3", s=5,alpha = .9)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_rna_dna_color_maxUMIrna.pdf')


# # the following to analyses create scatter that comapres rna counts between replicates and dna counts between replicates - but these counts should be normalized to barcode counts (see next analyses) (not main analyses)
# # correlation of rna cpm between replicates
# print("correlation of RNA cpm between replicates")
# for rep in [["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]:
    # print(rep)
    # corr_active_all = data_df.loc[(data_df["activity"] == "active"), f"RNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df["activity"] == "active"), f"RNA_pseudo_{rep[1]}"])
    # corr_active_5 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[1]}"])
    # corr_active_10 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"RNA_pseudo_{rep[1]}"])

    # print("all bcs:", corr_active_all)
    # print(">=5bcs:", corr_active_5)
    # print(">=10bcs:", corr_active_10)
    
    
    # values = np.vstack([data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"RNA_pseudo_{rep[0]}"], data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"RNA_pseudo_{rep[1]}"]])
    # kernel = stats.gaussian_kde(values)(values)
    # plt.clf()
    # sns.scatterplot(data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"RNA_pseudo_{rep[0]}", y = f"RNA_pseudo_{rep[1]}", c=kernel, s=10)
    # plt.title(f'{cells}, {library}{adaptor}, correlation {rep[0]} and {rep[1]} RNA cpm - active, >=5bcs')
    # plt.xlabel(f'RNA cpm {rep[0]}')
    # plt.ylabel(f'RNA cpm {rep[1]}')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_{rep[0]}_{rep[1]}_RNA.pdf')

# # correlation of dna cpm between replicates
# print("correlation of DNA cpm between replicates")
# for rep in [["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]:
    # print(rep)
    # corr_active_all = data_df.loc[(data_df["activity"] == "active"), f"DNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df["activity"] == "active"), f"DNA_pseudo_{rep[1]}"])
    # corr_active_5 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[1]}"])
    # corr_active_10 = data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[0]}"].corr(data_df.loc[(data_df[f"count_{rep[0]}"] >= 10)& (data_df[f"count_{rep[1]}"] >= 10)&(data_df["activity"] == "active"), f"DNA_pseudo_{rep[1]}"])

    # print("all bcs:", corr_active_all)
    # print(">=5bcs:", corr_active_5)
    # print(">=10bcs:", corr_active_10)
    
    
    # values = np.vstack([data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"DNA_pseudo_{rep[0]}"], data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"DNA_pseudo_{rep[1]}"]])
    # kernel = stats.gaussian_kde(values)(values)
    # plt.clf()
    # sns.scatterplot(data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"DNA_pseudo_{rep[0]}", y = f"DNA_pseudo_{rep[1]}", c=kernel, s=10)
    # plt.title(f'{cells}, {library}{adaptor}, correlation {rep[0]} and {rep[1]} DNA cpm - active, >=5bcs')
    # plt.xlabel(f'DNA cpm {rep[0]}')
    # plt.ylabel(f'DNA cpm {rep[1]}')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_{rep[0]}_{rep[1]}_DNA.pdf')



# # correlation between dna-dna and rna-rna  normalized to number of barcodes (maybe main analysis?)
# for molecule in ["DNA","RNA"]:
    # print(f"correlation of {molecule} cpm between replicates")
    # for rep in [["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]:
        # print(rep)
        
        # filtered_df = data_df[(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5) & (data_df["activity"] == "active")]
        
        # ratio1 = filtered_df[f"{molecule}_pseudo_{rep[0]}"]/filtered_df[f"count_{rep[0]}"]
        # ratio2 = filtered_df[f"{molecule}_pseudo_{rep[1]}"]/filtered_df[f"count_{rep[1]}"]
        
        # correlation = ratio1.corr(ratio2)
        # print(correlation)
        
        # values = np.vstack([ratio1, ratio2])
        # kernel = stats.gaussian_kde(values)(values)
        # plt.clf()
        # ax = sns.scatterplot(x = ratio1, y = ratio2, c=kernel, s=10)
        # ax.set_title(f'{cells}, {library}{adaptor}, correlation {rep[0]} and {rep[1]} {molecule} cpm - active, >=5bcs, normalized to #bc')
        # ax.set_xlabel(f'{molecule} cpm {rep[0]} / #bc')
        # ax.set_ylabel(f'{molecule} cpm {rep[1]} / #bc')
        # min_xy, max_xy = get_fig_xlim_ylim(ratio1,ratio2)
        # plt.xlim(right = max_xy)
        # plt.ylim(top = max_xy)
        # ax.set_aspect('equal', 'box')
        # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5bcs_{rep[0]}_{rep[1]}_{molecule}_bc_normalized.pdf')


# # correlation between possible barcodes and barcodes in RNA/DNA (not main analysis)
# print("correlation between possible barcodes and barcodes in RNA/DNA")

# for rep in ["rep1","rep2","rep3"]:
    # print(rep)
    # no_na = data_df.dropna(subset = [f'count_{rep}'])
    # corr_bcs = no_na['possible_bcs'].corr(no_na[f'count_{rep}'])
    # print(corr_bcs)
    
    # values = np.vstack([no_na[f"possible_bcs"], no_na[f"count_{rep}"]])
    # kernel = stats.gaussian_kde(values)(values)
    # plt.clf()
    # sns.scatterplot(no_na, x = f"possible_bcs", y = f"count_{rep}", c=kernel, s=10)
    # plt.title(f'{cells}, {library}{adaptor}, {rep}, correlation between #bcs in association and #bcs in rna/dna')
    # plt.xlabel(f'bcs association')
    # plt.ylabel(f'bcs rna/dna, {rep}')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_barcodes_Association_vs_rnadna_{rep}.pdf')
    # plt.clf()
    # sns.scatterplot(no_na, x = f"possible_bcs", y = f"count_{rep}", c=kernel, s=10)
    # plt.title(f'{cells}, {library}{adaptor}, {rep}, correlation between #bcs in association and #bcs in rna/dna, lim=200')
    # plt.ylim(0, 200)
    # plt.xlim(0, 200)
    # plt.xlabel(f'bcs association')
    # plt.ylabel(f'bcs rna/dna, {rep}')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_barcodes_Association_vs_rnadna_{rep}_lim200.pdf')

# sns.reset_defaults()

# histogram of RNA/DNA colored according to activity and stacked (MAIN ANALYSIS)
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True)
for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    sns.histplot(ax=axes[n], data=data_df, x=f"ratio_log_{rep}", hue = "activity", hue_order = ["non_active", "active"], multiple="stack", palette=["darkgrey", "darkorange"], edgecolor='none', alpha=1,binwidth=0.03 )
    # plt.xlim(-4, 4)
    # plt.ylim(0, 3200)
plt.suptitle("Distribution log RNA/DNA")
plt.savefig(f'{activity_folder_path}/RNA_DNA_distribution_stacked_{cells}_{library}{adaptor}_color.png', transparent=True, dpi=330)
plt.savefig(f'{activity_folder_path}/RNA_DNA_distribution_stacked_{cells}_{library}{adaptor}_color.pdf')

# the follwing analyses create histograms of rna and dna counts - have to think about it a bit more - should maybe be done on cpm nomralized counts insteas of raw counts (MAIN ANALYSIS)
# histogram of RNA
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True)
for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    sns.histplot(ax=axes[n], data=data_df, x=f"RNA_{rep}", edgecolor='none')
plt.suptitle("Distribution RNA")
plt.savefig(f'{activity_folder_path}/RNA_distribution_{cells}_{library}{adaptor}.png', transparent=True, dpi=330)
plt.savefig(f'{activity_folder_path}/RNA_distribution_{cells}_{library}{adaptor}.pdf')

# histogram of RNA - clipped to 5000
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True)
for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    curr_df = data_df.copy()
    curr_df[f"RNA_{rep}"] = curr_df[f"RNA_{rep}"].clip(upper=5000)
    sns.histplot(ax=axes[n], data=curr_df, x=f"RNA_{rep}", edgecolor='none')
plt.suptitle("Distribution RNA - clipped to 5000")
plt.savefig(f'{activity_folder_path}/RNA_distribution_{cells}_{library}{adaptor}_clipped5000.png', transparent=True, dpi=330)
plt.savefig(f'{activity_folder_path}/RNA_distribution_{cells}_{library}{adaptor}_clipped5000.pdf')

# histogram of DNA
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True)
for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    sns.histplot(ax=axes[n], data=data_df, x=f"DNA_{rep}", edgecolor='none')
plt.suptitle("Distribution DNA")
plt.savefig(f'{activity_folder_path}/DNA_distribution_{cells}_{library}{adaptor}.png', transparent=True, dpi=330)
plt.savefig(f'{activity_folder_path}/DNA_distribution_{cells}_{library}{adaptor}.pdf')

# distribution of alpha (MAIN ANALYSIS)
plt.clf()
sns.histplot(data=data_df, x="alpha", hue = "activity", edgecolor='none')
plt.title('Distribution of alpha')
plt.savefig(f'{activity_folder_path}/alpha_distribution_{cells}_{library}{adaptor}.pdf')
plt.savefig(f'{activity_folder_path}/alpha_distribution_{cells}_{library}{adaptor}.png', dpi=330)

plt.clf()
sns.histplot(data=data_df, x="alpha", hue = "activity", edgecolor='none')
plt.title('Distribution of alpha - xlim=20')
plt.xlim(0, 20)
plt.savefig(f'{activity_folder_path}/alpha_distribution_{cells}_{library}{adaptor}_xlim_20.pdf')
plt.savefig(f'{activity_folder_path}/alpha_distribution_{cells}_{library}{adaptor}_xlim_20.png', dpi=330)

# distribution of rna/dna ratio of scrambled (MAIN ANALYSIS)
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True)
for n, rep in enumerate(["rep1", "rep2", "rep3"]):
    data_scrambled = data_df[data_df["control_type"]== "scrambled_control"]
    sns.histplot(ax=axes[n], data=data_scrambled, x=f"ratio_log_{rep}")
plt.suptitle("Distribution log RNA/DNA scrambled sequences")
plt.savefig(f'{activity_folder_path}/scrambled_distribution_{cells}_{library}{adaptor}.pdf')
plt.savefig(f'{activity_folder_path}/scrambled_distribution_{cells}_{library}{adaptor}.png', dpi=330)

# get oligo coverage for each replicate (MAIN ANALYSIS)
print("oligo coverage")
for rep in ["rep1", "rep2", "rep3", "rep_comb"]:
    for bc in [1, 5, 10, 20]:
        count = (data_df[f"DNA_{rep}"] >= bc).sum()
        ratio = count/79812
        print(f"{rep}: count of oligos with >= {bc} DNA: {count}")
        print(f"{rep}: ratio of oligos with >= {bc} DNA: {ratio}")
