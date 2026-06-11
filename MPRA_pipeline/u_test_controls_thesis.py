# this script conducts a mann-whitney u test comapring the different controls groups
# in chondrocytes (L1a2) it compares all negative controls to all positive controls
# in neurons (l2a2) it compares ???

# cells=neurons
# library=L2
# adaptor=a2
# module load python/3.7.9
# bsub -q molgen-q -R "rusage[mem=4000]" python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/u_test_controls_thesis.py $library $adaptor $cells

import pandas as pd
from scipy import stats
from scipy.stats import mannwhitneyu
import sys
import numpy as np

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]


activity_folder_path = f'./{cells}/{library}{adaptor}/output/activity_after_filter'
comb_df_path = f'{activity_folder_path}/comb_df_adjusted_fdr.csv'
data_df = pd.read_csv(comb_df_path)

pos = data_df[data_df['control_type'].str.contains('NegCtrl_active_not_diff|PosCtrl_chondrocyte_active|PosCtrl_diff|PosCtrl_neuron_active|PosCtrl_osteoblast_active', case=False, regex=True)]['alpha']
print(np.nanmean(pos))

neg = data_df[data_df['control_type'].str.contains('NegCtrl_not_active|scrambled_control|NegCtrl_non_SCREEN', case=False, regex=True)]['alpha']
print(np.nanmean(neg))

U1, p = mannwhitneyu(pos, neg, alternative = "greater")

print(p)