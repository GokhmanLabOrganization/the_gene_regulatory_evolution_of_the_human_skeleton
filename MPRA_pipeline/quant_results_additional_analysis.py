# this script does to additional analysis after the comb_df is created
# (we thought of them specifically to investigate the issue with the neurons, but obviously they can also be done on the chondrocytes)
# analysis 1: correlation of sequences that are positive controls for activity
# analysis 2: correlation between RNA/DNA and alpha
# analysis 3: calculate correlations between the replicates with higher DNA count thresholds
# analysis 4: testing correlatiom between replicates when switching rna and dna counts
# analysis 5: testing different dna thresholds - plotting correlation vs coverage
# analysis 6: histogram of dna counts per oligo (dna rep comb) and reverse cumulative density plot

import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from scipy import stats
from scipy.stats import pearsonr

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

activity_folder_path = f'./{cells}/{library}{adaptor}/output/activity_after_filter'
comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv'

data_df = pd.read_csv(comb_df_path)

# # analysis 1: correlation of sequences that are positive controls for activity

# for control_type in ["PosCtrl_osteoblast_active", "PosCtrl_neuron_active"]:
    # print(control_type)
    # plt.clf()
    # fig, axes = plt.subplots(1, 3, subplot_kw=dict(box_aspect=1))

    # fig.suptitle(f'Correlation of {control_type} between replicates')
    # for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
        # sns.scatterplot(ax=axes[n], data = data_df[(data_df[f"control_type"] == f"{control_type}")&(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5)], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", s=8)
        # print(rep)
        # corr_active_all = data_df.loc[(data_df[f"control_type"] == f"{control_type}"), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"control_type"] == f"{control_type}"), f"ratio_log_{rep[1]}"])
        # corr_active_5 = data_df.loc[(data_df[f"control_type"] == f"{control_type}")&(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5), f"ratio_log_{rep[0]}"].corr(data_df.loc[(data_df[f"control_type"] == f"{control_type}")&(data_df[f"count_{rep[0]}"] >= 5)& (data_df[f"count_{rep[1]}"] >= 5), f"ratio_log_{rep[1]}"])
        # print("all bcs:", corr_active_all)
        # print(">=5bcs:", corr_active_5)
        # # axes[n].axis('equal')
    # plt.tight_layout()
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_5bcs_{control_type}.pdf')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_5bcs_{control_type}.png', dpi=330)
    
# # analysis 2: correlation between RNA/DNA and alpha
# data_df_no_na = data_df.dropna(subset=['alpha', 'ratio_log_rep_comb'], how='any')
# print(data_df.shape[0])
# print(data_df_no_na.shape[0])
# data_df_no_na['ratio_log_rep_comb_no_log2'] = 2 ** data_df_no_na["ratio_log_rep_comb"]
# values12 = np.vstack([data_df_no_na["alpha"], data_df_no_na["ratio_log_rep_comb_no_log2"]])
# kernel12 = stats.gaussian_kde(values12)(values12)
# plt.clf()
# plt.figure(figsize=(7,7))
# ax = sns.scatterplot(data_df_no_na, x = "alpha", y = "ratio_log_rep_comb_no_log2", c=kernel12, s=10)
# plt.title(f'{cells}, {library}{adaptor}, alpha vs RNA/DNA')
# plt.xlabel('alpha')
# plt.ylabel('RNA/DNA')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_alpha_ratio.pdf')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_alpha_ratio.png', dpi=330)


# correlation = data_df_no_na['alpha'].corr(data_df_no_na["ratio_log_rep_comb_no_log2"])
# print("corr alpha, rna/dna:", correlation)

# analysis 3: calculate correlations between the replicates with higher DNA count thresholds

data_df_no_na_ratio_rep1_rep2 = data_df.dropna(subset=['ratio_log_rep1', 'ratio_log_rep2'], how='any')
print(data_df.shape[0])
print(data_df_no_na_ratio_rep1_rep2.shape[0])

for activity in ['all', 'active']:
    print(activity)
    if activity == "all":
        data_df_no_na_ratio_rep1_rep2 = data_df_no_na_ratio_rep1_rep2.copy()
    elif activity == "active":
        data_df_no_na_ratio_rep1_rep2 = data_df_no_na_ratio_rep1_rep2[data_df_no_na_ratio_rep1_rep2["activity"] == 'active'].copy()
        
    # plt.clf()
    # fig, axes = plt.subplots(1,4, sharey=True, sharex=True, figsize=[10, 8])
    for n, DNA_thres in enumerate([5, 500, 1000, 2000]):
        print(DNA_thres)
        data_df_sub = data_df_no_na_ratio_rep1_rep2[data_df_no_na_ratio_rep1_rep2["DNA_rep_comb"] > DNA_thres].copy()
        print(data_df_sub.shape[0])

        corr = data_df_sub.loc[:, f"ratio_log_rep1"].corr(data_df_sub.loc[:, f"ratio_log_rep2"])
        print(corr)
        print(pearsonr(data_df_sub.loc[:, f"ratio_log_rep1"], data_df_sub.loc[:, f"ratio_log_rep2"]))
        
        # values = np.vstack([data_df_sub[f"ratio_log_rep1"], data_df_sub[f"ratio_log_rep2"]])
        # kernel = stats.gaussian_kde(values)(values)
        # sns.scatterplot(ax=axes[n], data = data_df_sub, x = f"ratio_log_rep1", y = f"ratio_log_rep2", c=kernel, s=1, linewidth=0)
        # axes[n].set(xlabel=f'rep1', ylabel=f'rep2')
        # # axes[n].set_box_aspect(1)
        # axes[n].set_aspect('equal')
        # axes[n].set_title(f"DNA threshold: {DNA_thres}")
        # # axes[n].yaxis.get_label().set_visible(True)
    # plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_rna_dna_high_DNA_thresholds_{activity}.pdf')
    # plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_rna_dna_high_DNA_thresholds_{activity}.png', dpi=330)
    
# # analysis 4: testing correlatiom between replicates when switching rna and dna counts

# cpm = pd.read_csv(f'./{cells}/{library}{adaptor}/output/activity_after_filter/ratio_after_filter.csv')


# cpm["ratio_switched_rep1"] = np.log2(cpm["DNA_cpm_rep1"]/cpm["RNA_cpm_rep1"])
# cpm["ratio_switched_rep2"] = np.log2(cpm["DNA_cpm_rep2"]/cpm["RNA_cpm_rep2"])

# cpm_no_na = cpm.dropna(subset=['ratio_switched_rep1', 'ratio_switched_rep2'], how='any')
# print(cpm.shape[0])
# print(cpm_no_na.shape[0])

# plt.clf()

# values = np.vstack([cpm_no_na[f"ratio_switched_rep1"], cpm_no_na[f"ratio_switched_rep2"]])
# kernel = stats.gaussian_kde(values)(values)
# sns.scatterplot(data = cpm_no_na, x = f"ratio_switched_rep1", y = f"ratio_switched_rep2", c=kernel, s=1, linewidth=0)
# plt.axes().set(xlabel=f'rep1', ylabel=f'rep2')
# plt.axes().set_aspect('equal')
# plt.axes().set_title('RNA DNA switched')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_corr_ratio_rna_dna_switched.pdf')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_corr_ratio_rna_dna_switched.png', dpi=330)

# # analysis 5: testing different dna thresholds - plotting correlation vs coverage
# correlations_all = []
# correlations_active = []
# coverage_all = []
# coverage_active = []


# for DNA_thres in [5, 100, 200, 300, 400, 500, 600, 1000]:
    # data_thres = data_df[data_df["DNA_rep_comb"] > DNA_thres].copy()
    # corr_all = data_thres.loc[:, f"ratio_log_rep1"].corr(data_thres.loc[:, f"ratio_log_rep2"])
    # corr_active = data_thres.loc[(data_df["activity"] == "active"), f"ratio_log_rep1"].corr(data_thres.loc[(data_df["activity"] == "active"), f"ratio_log_rep2"])
    
    # cov_all = data_thres.shape[0]/data_df.shape[0]
    # cov_active = data_thres[data_thres["activity"] == "active"].shape[0]/data_df[data_df["activity"] == "active"].shape[0]
    
    # correlations_all.append(corr_all)
    # correlations_active.append(corr_active)
    # coverage_all.append(cov_all)
    # coverage_active.append(cov_active)
    
# print(correlations_all)
# print(correlations_active)
# print(coverage_all)
# print(coverage_active)

# plt.clf()
# fig, axs = plt.subplots(1, 3, figsize=[15, 8])
# axs[0].plot(correlations_all, coverage_all, 'o-')
# axs[1].plot(correlations_active, coverage_active, 'o-')
# axs[2].plot(correlations_active, coverage_all, 'o-')
# axs[0].set_xlim([0,1])
# axs[1].set_xlim([0,1])
# axs[2].set_xlim([0,1])
# axs[0].set_ylim([0,1])
# axs[1].set_ylim([0,1])
# axs[2].set_ylim([0,1])
# axs[0].set_xlabel("correlation_all")
# axs[1].set_xlabel("correlation_active")
# axs[2].set_xlabel("correlation_active")
# axs[0].set_ylabel("coverage_all")
# axs[1].set_ylabel("coverage_active")
# axs[2].set_ylabel("coverage_all")

# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_correlation_coverage_dna_threshold.pdf')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_correlation_coverage_dna_threshold.png', dpi=330)

# # analysis 6: histogram of dna counts per oligo (dna rep comb) and reverse cumulative density plot
# plt.clf()
# curr_df = data_df.copy()
# curr_df[f"DNA_rep_comb"] = curr_df[f"DNA_rep_comb"].clip(upper=7500)
# sns.histplot(data=curr_df, x=f"DNA_rep_comb", edgecolor='none')
# plt.suptitle("Distribution DNA - clipped to 7500")
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_DNA_distribution_{library}{adaptor}_clipped7500.png', transparent=True, dpi=330)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_DNA_distribution_{library}{adaptor}_clipped7500.pdf')

# plt.clf()
# curr_df = data_df.copy()
# sns.ecdfplot(data=curr_df, x=f"DNA_rep_comb", complementary=True)
# plt.suptitle("reverse cdf DNA")
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_reverse_cdf_DNA_{library}{adaptor}.png', transparent=True, dpi=330)
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_reverse_cdf_DNA_{library}{adaptor}.pdf')