import pandas as pd
import seaborn as sns
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import sys
import numpy as np


cells = sys.argv[1]

list_ctrl_df = []

for library in ["L1","L2","L3","L4"]:
    for adaptor in ["a1","a2","a3"]:
        if library == 'L4' and adaptor == 'a3':
            continue
        activity_folder_path = f'./{cells}/{library}{adaptor}/output/activity_after_filter'
        comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv' # NM 30-10-2024 Changed com_df.csv to comb_df_adjusted_fdr.csv 
        data_df = pd.read_csv(comb_df_path, usecols=['oligo', 'alpha', 'fdr.mad', 'ratio_log_rep_comb', 'control_type', 'activity'])
        data_df_ctrls = data_df[(data_df["control_type"] != "No control")]
        #data_df_ctrls = data_df[data_df["control_type"].str.contains("Ctrl")]

        data_df_ctrls["oligo"] = data_df_ctrls["oligo"].str[:-6]
        list_ctrl_df.append(data_df_ctrls)

print(list_ctrl_df[0].head().to_string())

merged_ctrls = list_ctrl_df[0].merge(list_ctrl_df[1],on='oligo',how="outer",suffixes=("_L1a1",None)).merge(
list_ctrl_df[2],on='oligo',how="outer",suffixes=("_L1a2",None)).merge(
list_ctrl_df[3],on="oligo",how="outer",suffixes=("_L1a3",None)).merge(
list_ctrl_df[4],on='oligo',how="outer",suffixes=("_L2a1",None)).merge(
list_ctrl_df[5],on='oligo',how="outer",suffixes=("_L2a2",None)).merge(
list_ctrl_df[6],on="oligo",how="outer",suffixes=("_L2a3",None)).merge(
list_ctrl_df[7],on='oligo',how="outer",suffixes=("_L3a1",None)).merge(
list_ctrl_df[8],on="oligo",how="outer",suffixes=("_L3a2",None)).merge(
list_ctrl_df[9],on="oligo",how="outer",suffixes=("_L3a3",None)).merge(
list_ctrl_df[10],on="oligo",how="outer",suffixes=("_L4a1","_L4a2"))

merged_ctrls.to_csv(f'./additional/analyze_controls/{cells}/controls_df.csv', header=True, index = False)


print(merged_ctrls.head().to_string())

alpha = merged_ctrls.filter(regex='alpha')

corr_alpha = alpha.corr()
plt.clf()
sns.heatmap(corr_alpha, cmap="Blues", annot=True)
plt.tight_layout()
plt.savefig(f'./additional/analyze_controls/{cells}/alpha_corr_heatmap.pdf')
plt.savefig(f'./additional/analyze_controls/{cells}/alpha_corr_heatmap.png', dpi=330)


# # add the neurons to the heatmap - this is adjusted only to the three neurons batchrs that I have at the moment - can do this differently when we have all batches #KL 21.09.23
# list_neurons = []
# for library_adaptor in ["L1a1", "L1a2", "L3a3"]:
    # activity_folder_path = f'./neurons/{library_adaptor}/output/activity_after_filter'
    # comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv' # NM 30-10-2024 Changed com_df.csv to comb_df_adjusted_fdr.csv 
    # data_df = pd.read_csv(comb_df_path, usecols=['oligo', 'alpha', 'fdr.mad', 'ratio_log_rep_comb', 'control_type', 'activity'])
    # data_df_ctrls = data_df[(data_df["control_type"] != "No control")]
    # data_df_ctrls["oligo"] = data_df_ctrls["oligo"].str[:-6]
    # list_neurons.append(data_df_ctrls)
    
# merged_neurons = list_neurons[0].merge(list_neurons[1],on='oligo',how="outer",suffixes=("_L1a1",None)).merge(list_neurons[2],on='oligo',how="outer",suffixes=("_L1a2","_L3a3"))
# merged_neurons_chondro = merged_ctrls.merge(merged_neurons,on='oligo',how="outer",suffixes=(None,"_neurons"))
# alpha_neurons_chondro = merged_neurons_chondro.filter(regex='alpha')

# corr_alpha_neurons_chondro = alpha_neurons_chondro.corr()
# plt.clf()
# sns.heatmap(corr_alpha_neurons_chondro, cmap="Blues", annot=True, annot_kws={"size":4})
# plt.tight_layout()
# plt.savefig(f'./additional/analyze_controls/alpha_corr_heatmap_chondro_neurons.pdf')
# plt.savefig(f'./additional/analyze_controls/alpha_corr_heatmap_chondro_neurons.png', dpi=330)

# plt.clf()
# fig, axis = plt.subplots(9,1, sharex=True,figsize=(10, 18))
# num_bins=50
# mini=np.nanmin(alpha.to_numpy())
# maxi=np.nanmax(alpha.to_numpy())
# print(mini, maxi)
# alpha.hist(ax=axis, bins=num_bins, range=(mini,maxi))
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_hist.pdf')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_hist.png', dpi=330)

# plt.clf()
# fig, axis = plt.subplots(9,1, sharex=True,figsize=(10, 18))
# num_bins=100
# mini=np.nanmin(alpha.to_numpy())
# maxi=np.nanmax(alpha.to_numpy())
# alpha.hist(ax=axis, bins=num_bins, range=(mini,maxi))
# plt.xlim(0, 20)
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_hist_xlim20.pdf')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_hist_xlim20.png', dpi=330)

# plt.clf()
# fig, axes = plt.subplots(9, 9, sharex=True, sharey=True,figsize=(10, 10))
# for n in range(9):
    # for m in range(9):
        # sns.scatterplot(ax=axes[n, m], data = alpha, x = alpha.columns[m], y = alpha.columns[n], s=5, alpha = .9)
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter.pdf')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter.png', dpi=330)

# plt.clf()
# fig, axes = plt.subplots(9, 9, sharex=True, sharey=True,figsize=(10, 10))
# for n in range(9):
    # for m in range(9):
        # sns.scatterplot(ax=axes[n, m], data = alpha, x = alpha.columns[m], y = alpha.columns[n], s=5, alpha = .9)
# plt.xscale('log')
# plt.yscale('log')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter_log.pdf')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter_log.png', dpi=330)

# plt.clf()
# sns.scatterplot(data = alpha, x = "alpha_L1a1", y = "alpha_L1a2", s=5, alpha = .9)
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter_L1a1_L1a2.pdf')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter_L1a1_L1a2.png', dpi=330)

# plt.clf()
# sns.scatterplot(data = alpha, x = "alpha_L3a2", y = "alpha_L1a2", s=5, alpha = .9)
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter_L3a2_L1a2.pdf')
# plt.savefig(f'./additional/analyze_controls/{cells}/alpha_scatter_L3a2_L1a2.png', dpi=330)

# print("all")
# for library in ["L1","L2","L3"]:
    # for adaptor in ["a1","a2","a3"]:
        # num_active = (merged_ctrls[f"activity_{library}{adaptor}"] == "active").sum()
        # num_not_active = (merged_ctrls[f"activity_{library}{adaptor}"] == "non_active").sum()
        # print(f"{library}{adaptor}")
        # print(f'# active: {num_active}')
        # print(f'# not active: {num_not_active}')
        # print(f'% active: {num_active/(num_active+num_not_active)}')
        
# print("no_na")
# merged_ctrls_no_na = merged_ctrls.dropna()
# for library in ["L1","L2","L3"]:
    # for adaptor in ["a1","a2","a3"]:
        # num_active = (merged_ctrls_no_na[f"activity_{library}{adaptor}"] == "active").sum()
        # num_not_active = (merged_ctrls_no_na[f"activity_{library}{adaptor}"] == "non_active").sum()
        # print(f"{library}{adaptor}")
        # print(f'# active: {num_active}')
        # print(f'# not active: {num_not_active}')
        # print(f'% active: {num_active/(num_active+num_not_active)}')