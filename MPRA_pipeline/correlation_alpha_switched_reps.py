import pandas as pd
import sys
import numpy as np
from scipy import stats
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

comb_df_path_switch = f'./additional/switching_replicates/{cells}/{library}{adaptor}/output/activity/comb_df.csv'
comb_df_path = f'./{cells}/{library}{adaptor}/output/activity/comb_df.csv'

switched = pd.read_csv(comb_df_path_switch)
not_switched = pd.read_csv(comb_df_path)

merged = not_switched.merge(switched, on='oligo', suffixes=('_not_switched', '_switched'))

corr = merged['alpha_not_switched'].corr(merged['alpha_switched'])
print(corr)

no_na = merged.dropna(subset = ['alpha_not_switched', 'alpha_switched'])

corr = no_na['alpha_not_switched'].corr(no_na['alpha_switched'])
print(corr)

values = np.vstack([no_na[f"alpha_not_switched"], no_na[f"alpha_switched"]])
kernel = stats.gaussian_kde(values)(values)
plt.clf()
sns.scatterplot(no_na, x = f"alpha_not_switched", y = f"alpha_switched", c=kernel, s=10)
plt.title(f'{cells}, {library}{adaptor}, correlation between alphas, log')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(f'alpha - reps not switched')
plt.ylabel(f'alpha - reps not switched')
plt.savefig(f'./additional/switching_replicates/{cells}/{library}{adaptor}/output/activity/{cells}_{library}{adaptor}_scatter_alpha_switched_reps_log.pdf')
