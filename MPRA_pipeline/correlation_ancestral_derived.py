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

data_ancestral = data_df[data_df['oligo'].str.contains("ancestral")]
print(data_ancestral.shape[0])
data_derived = data_df[data_df['oligo'].str.contains("derived")]
print(data_derived.shape[0])

data_ancestral['oligo_pure'] = data_ancestral['oligo'].str.replace("_ancestral", '')
data_derived['oligo_pure'] = data_derived['oligo'].str.replace("_derived", '')


data_merged = data_ancestral.merge(data_derived, how="inner", on = "oligo_pure", suffixes=('_ancestral', '_derived'))

# plt.clf()
# fig, axes = plt.subplots(1,5, figsize=(20, 4))
# for n, rep in enumerate(["rep1","rep2","rep3", "rep_comb", "alpha"]):
    # if rep == "alpha":
        # print(rep)
        # corr = data_merged[f'{rep}_ancestral'].corr(data_merged[f'{rep}_derived'])
        # print(corr)
        
        # # values = np.vstack([data_merged[f'{rep}_ancestral'], data_merged[f'{rep}_derived']])
        # # kernel = stats.gaussian_kde(values)(values)

        # sns.scatterplot(ax=axes[n], data=data_merged, x = f"{rep}_ancestral", y = f"{rep}_derived", s=10)
        # axes[n].set_aspect('equal')
    # else:
        # print(rep)
        # corr = data_merged[f'ratio_log_{rep}_ancestral'].corr(data_merged[f'ratio_log_{rep}_derived'])
        # print(corr)
        # data_merged_temp = data_merged.dropna(subset = [f'ratio_log_{rep}_ancestral',f'ratio_log_{rep}_derived'])
        # print(pearsonr(data_merged_temp[f'ratio_log_{rep}_ancestral'], data_merged_temp[f'ratio_log_{rep}_derived']))
        # values = np.vstack([data_merged[f'ratio_log_{rep}_ancestral'], data_merged[f'ratio_log_{rep}_derived']])
        # kernel = stats.gaussian_kde(values)(values)

        # sns.scatterplot(ax=axes[n], data=data_merged, x = f"ratio_log_{rep}_ancestral", y = f"ratio_log_{rep}_derived", s=10)
        # axes[n].set_aspect('equal')
        
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation between ancestral and derived allele')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral.pdf')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral.png', dpi=330)


# print("active sequences")
# data_merged_active = data_merged[(data_merged['activity_adjusted_ancestral']== "active") | (data_merged['activity_adjusted_derived'] == "active")]

# plt.clf()
# fig, axes = plt.subplots(1,5, figsize=(20, 4))
# for n, rep in enumerate(["rep1","rep2","rep3", "rep_comb", "alpha"]):
    # if rep == "alpha":
        # print(rep)
        # corr = data_merged_active[f'{rep}_ancestral'].corr(data_merged_active[f'{rep}_derived'])
        # print(corr)
        
        # values = np.vstack([data_merged_active[f'{rep}_ancestral'], data_merged_active[f'{rep}_derived']])
        # kernel = stats.gaussian_kde(values)(values)

        # sns.scatterplot(ax=axes[n], data=data_merged_active, x = f"{rep}_ancestral", y = f"{rep}_derived", s=10)
        # axes[n].set_aspect('equal')
    # else:
        # print(rep)
        # corr = data_merged_active[f'ratio_log_{rep}_ancestral'].corr(data_merged_active[f'ratio_log_{rep}_derived'])
        # print(corr)
        # data_merged_active_temp = data_merged_active.dropna(subset = [f'ratio_log_{rep}_ancestral',f'ratio_log_{rep}_derived'])
        # print(pearsonr(data_merged_active_temp[f'ratio_log_{rep}_ancestral'], data_merged_active_temp[f'ratio_log_{rep}_derived']))
        
        # values = np.vstack([data_merged_active[f'ratio_log_{rep}_ancestral'], data_merged_active[f'ratio_log_{rep}_derived']])
        # kernel = stats.gaussian_kde(values)(values)

        # sns.scatterplot(ax=axes[n], data=data_merged_active, x = f"ratio_log_{rep}_ancestral", y = f"ratio_log_{rep}_derived", s=10)
        # axes[n].set_aspect('equal')
        
# plt.suptitle(f'{cells}, {library}{adaptor}, correlation between ancestral and derived allele for active oligos')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral_active.pdf')
# plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral_active.png', dpi=330)

# for thesis: 

plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True, figsize=(10, 4))
for n, rep in enumerate(["rep1","rep2","rep3"]):
    print(rep)
    corr = data_merged[f'ratio_log_{rep}_ancestral'].corr(data_merged[f'ratio_log_{rep}_derived'])
    print(corr)
    data_merged_temp = data_merged.dropna(subset = [f'ratio_log_{rep}_ancestral',f'ratio_log_{rep}_derived'])
    print(pearsonr(data_merged_temp[f'ratio_log_{rep}_ancestral'], data_merged_temp[f'ratio_log_{rep}_derived']))
    values = np.vstack([data_merged_temp[f'ratio_log_{rep}_ancestral'], data_merged_temp[f'ratio_log_{rep}_derived']])
    kernel = stats.gaussian_kde(values)(values)

    sns.scatterplot(ax=axes[n], data=data_merged_temp, x = f"ratio_log_{rep}_ancestral", y = f"ratio_log_{rep}_derived", c=kernel, s=1, linewidth=0)
    axes[n].set_aspect('equal')
        
plt.suptitle(f'{cells}, {library}{adaptor}, correlation between ancestral and derived allele')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral_for_thesis.pdf')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral_for_thesis.png', dpi=330)


print("active sequences")
data_merged_active = data_merged[(data_merged['activity_adjusted_ancestral']== "active") | (data_merged['activity_adjusted_derived'] == "active")]

plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True, figsize=(10, 4))
for n, rep in enumerate(["rep1","rep2","rep3"]):
    print(rep)
    corr = data_merged_active[f'ratio_log_{rep}_ancestral'].corr(data_merged_active[f'ratio_log_{rep}_derived'])
    print(corr)
    data_merged_active_temp = data_merged_active.dropna(subset = [f'ratio_log_{rep}_ancestral',f'ratio_log_{rep}_derived'])
    print(pearsonr(data_merged_active_temp[f'ratio_log_{rep}_ancestral'], data_merged_active_temp[f'ratio_log_{rep}_derived']))
    
    values = np.vstack([data_merged_active_temp[f'ratio_log_{rep}_ancestral'], data_merged_active_temp[f'ratio_log_{rep}_derived']])
    kernel = stats.gaussian_kde(values)(values)

    sns.scatterplot(ax=axes[n], data=data_merged_active_temp, x = f"ratio_log_{rep}_ancestral", y = f"ratio_log_{rep}_derived", c=kernel, s=1, linewidth=0)
    axes[n].set_aspect('equal')
        
plt.suptitle(f'{cells}, {library}{adaptor}, correlation between ancestral and derived allele for active oligos')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral_active_for_thesis.pdf')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_derived_ancestral_active_for_thesis.png', dpi=330)