# cells=neurons
# library=L2
# adaptor=a2
# module load python/3.7.9
# bsub -q molgen-q -R "rusage[mem=4000]" python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/scatter_for_thesis.py $library $adaptor $cells


import pandas as pd
import sys
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

activity_folder_path = f'./{cells}/{library}{adaptor}/output/activity_after_filter'
comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv'
data_df = pd.read_csv(comb_df_path)


# correlation between reps
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True, figsize=[10, 8])
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
    
    values = np.vstack([data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)][f"ratio_log_{rep[0]}"], data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)][f"ratio_log_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    sns.scatterplot(ax=axes[n], data = data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", c=kernel, s=1, linewidth=0)
    axes[n].set(xlabel=f'{rep[0]}', ylabel=f'{rep[1]}')
    axes[n].set_aspect('equal')
    axes[n].yaxis.get_label().set_visible(True)
plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - >=5bcs')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_5_DNA_rna_dna_for_thesis.pdf')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_5_DNA_rna_dna_for_thesis.png', dpi=330)

# correlation for only active sequences (here activity is defined according to the fdr that is done on all seqeunces - activity column and not activity_adjusted)
plt.clf()
fig, axes = plt.subplots(1,3, sharey=True, sharex=True, figsize=[10, 8])
for n, rep in enumerate([["rep1", "rep2"], ["rep2", "rep3"], ["rep3", "rep1"]]):
   
    values = np.vstack([data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5) & (data_df["activity"] == "active")][f"ratio_log_{rep[0]}"], data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)& (data_df["activity"] == "active")][f"ratio_log_{rep[1]}"]])
    kernel = stats.gaussian_kde(values)(values)
    sns.scatterplot(ax=axes[n], data = data_df[(data_df[f"DNA_{rep[0]}"] >= 5)& (data_df[f"DNA_{rep[1]}"] >= 5)& (data_df["activity"] == "active")], x = f"ratio_log_{rep[0]}", y = f"ratio_log_{rep[1]}", c=kernel, s=1, linewidth=0)
    axes[n].set(xlabel=f'{rep[0]}', ylabel=f'{rep[1]}')
    axes[n].set_aspect('equal')
    axes[n].yaxis.get_label().set_visible(True)
plt.suptitle(f'{cells}, {library}{adaptor}, correlation rna/dna, log2 - active, >=5bcs')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5_DNA_rna_dna_for_thesis.pdf')
plt.savefig(f'{activity_folder_path}/{cells}_{library}{adaptor}_scatter_active_5_DNA_rna_dna_for_thesis.png', dpi=330)