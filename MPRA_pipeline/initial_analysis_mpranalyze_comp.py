# creates histogram of acitve/differntial active fold change

import pandas as pd
import sys
import os
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns

cells = sys.argv[1]


output = pd.read_csv(f"/home/labs/davidgo/Collaboration/humanMPRA/{cells}/comparative_analysis_combined/mpranalyze_comp_res_fdr_not_complete.csv", header=0)


# histogram of RNA/DNA colored according to activity and stacked
plt.clf()
sns.displot(data=output, x="logFC", hue = "differential.full", multiple="stack", palette=["orange", "r"], edgecolor="none", alpha=1)
plt.savefig(f'./{cells}/comparative_analysis_combined/logFC_distribution_stacked.png', transparent=True, dpi=1200)
plt.savefig(f'./{cells}/comparative_analysis_combined/logFC_distribution_stacked.pdf')