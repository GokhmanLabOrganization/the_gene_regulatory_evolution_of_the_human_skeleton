# read in info_df, oligo fasta, chengyus fasta, chengyus results
# turn fastas into df
# remove adaptor from oligo fasta
# merge sequence on info df using oligo name
# merge chengyus name usign sequence
# merge chegyus results using chengyus name
# get numbers of overlapping sequences
# calculate correlation

import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# read info_df
info_df = pd.read_csv("/home/labs/davidgo/Collaboration/humanMPRA/additional/mpranalyze_quantitative/neurons_L1a1/info_df.csv")

# read chengyus results
c_results = pd.read_csv("/home/labs/davidgo/Collaboration/humanMPRA/additional/chengyu/normalized_ratio_pval.csv")

# read in oligo fasta and turn into df
oligo_dict = defaultdict(list)
with open('/home/labs/davidgo/Collaboration/humanMPRA/oligo_fasta/L1a1.fasta') as fp:
  for record in SeqIO.parse(fp,"fasta"):
    description = record.description
    sequence = str(record.seq)
    oligo_dict['description'].append(description)
    oligo_dict['sequence'].append(sequence)

oligo_fasta_df = pd.DataFrame.from_dict(oligo_dict)

print(oligo_fasta_df.head())

# remove last character from our sequence - sequence extension was done differently
oligo_fasta_df["sequence"] = oligo_fasta_df["sequence"].str[:-1]


# remove adaptor sequence from oligo_fasta_df - not needed: .seq already discards the adaptor apparently
# oligo_fasta_df["sequence"] = oligo_fasta_df["sequence"].str[15:-15]

# read in chengyus fasta
c_fasta_dict = defaultdict(list)
with open('/home/labs/davidgo/Collaboration/humanMPRA/additional/chengyu/lentiMPRA_constructs_270bp_nodup_regex_v210830_final_adapter+ctrl.txt') as fp:
  for record in SeqIO.parse(fp,"fasta"):
    c_name = record.description
    sequence = str(record.seq).lower()
    c_fasta_dict['c_name'].append(c_name)
    c_fasta_dict['sequence'].append(sequence)
    
c_fasta_df = pd.DataFrame.from_dict(c_fasta_dict)
# remove first character from changyus sequence - apparently sequence extension was performed differently
c_fasta_df["sequence"] = c_fasta_df["sequence"].str[1:]


print(c_fasta_df.head())

# merge sequence on info_df using oligo_fasta_df
info_df_seq = info_df.merge(oligo_fasta_df, left_on="oligo", right_on="description")
#print(info_df_seq.shape)

# merge c_name on info_df using c_fasta_df
info_df_seq_cname = info_df_seq.merge(c_fasta_df, on="sequence")
#print(info_df_seq_cname.head().to_string())

# merge c_results on info_df using c_name
info_df_seq_cname_cresults = info_df_seq_cname.merge(c_results, left_on = "c_name", right_on="insert")
#print(info_df_seq_cname_cresults.head().to_string())

# calculate log2 for chengyus ratio
info_df_seq_cname_cresults["mean_ratio_log2"] = np.log2(info_df_seq_cname_cresults['mean_ratio'])

# create column for top and bottom hue
info_df_seq_cname_cresults["top_bottom"] = info_df_seq_cname_cresults["c_name"].str[:24]

# save info_df_seq_cname_cresults
info_df_seq_cname_cresults.to_csv("/home/labs/davidgo/Collaboration/humanMPRA/additional/chengyu/comb_df.csv", index=False)

plt.clf()
sns.scatterplot(info_df_seq_cname_cresults, x = "ratio_log", y = "mean_ratio_log2", hue = "top_bottom")
plt.title(f'Correlation of expression in human MPRA and Chengyus experiment')
plt.xlabel('log2 RNA/DNA human MPRA')
plt.ylabel('log2 RNA/DNA Chengyu')
plt.savefig(f'/home/labs/davidgo/Collaboration/humanMPRA/additional/chengyu/neurons_correlation_chengyu_log2.pdf')

corr = info_df_seq_cname_cresults["ratio_log"].corr(info_df_seq_cname_cresults["mean_ratio"])
print(corr)
corr_new = info_df_seq_cname_cresults["ratio_log"].corr(info_df_seq_cname_cresults["mean_ratio_log2"])
print(corr_new)