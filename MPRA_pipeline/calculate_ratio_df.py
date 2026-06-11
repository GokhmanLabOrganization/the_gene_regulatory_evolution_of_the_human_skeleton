import pandas as pd
import numpy as np
import sys


library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

#Obtaining the read counts with soft coding (replacing the numbers above), Nadav 06-10-2024

UMI_all_reps = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t')
UMI_all_reps.set_index('oligo_bc', inplace=True)

rep1DNA=UMI_all_reps[f"DNA_{cells}_rep1"].sum()
rep1RNA=UMI_all_reps[f"RNA_{cells}_rep1"].sum()
rep2DNA=UMI_all_reps[f"DNA_{cells}_rep2"].sum()
rep2RNA=UMI_all_reps[f"RNA_{cells}_rep2"].sum()
rep3DNA=UMI_all_reps[f"DNA_{cells}_rep3"].sum()
rep3RNA=UMI_all_reps[f"RNA_{cells}_rep3"].sum()


# acess UMI counts file
rep_all = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t', header=0)

# access replicate UMI counts files
# rep1 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep1.txt', sep='\t', header=0, names=["oligo_bc", "RNA", "DNA"])
# rep2 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep2.txt', sep='\t', header=0, names=["oligo_bc", "RNA", "DNA"])
# rep3 = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_rep3.txt', sep='\t', header=0, names=["oligo_bc", "RNA", "DNA"])

#divide rep_all into reps
rep1 = rep_all[["oligo_bc", f"RNA_{cells}_rep1", f"DNA_{cells}_rep1"]]
rep1.rename(columns={f"RNA_{cells}_rep1": "RNA", f"DNA_{cells}_rep1":"DNA"}, inplace=True)
rep2 = rep_all[["oligo_bc", f"RNA_{cells}_rep2", f"DNA_{cells}_rep2"]]
rep2.rename(columns={f"RNA_{cells}_rep2": "RNA", f"DNA_{cells}_rep2":"DNA"}, inplace=True)
rep3 = rep_all[["oligo_bc", f"RNA_{cells}_rep3", f"DNA_{cells}_rep3"]]
rep3.rename(columns={f"RNA_{cells}_rep3": "RNA", f"DNA_{cells}_rep3":"DNA"}, inplace=True)

#sum counts across replicates
rep_all['RNA'] = rep_all[[col for col in rep_all.columns if 'RNA' in col]].sum(axis=1)
rep_all['DNA'] = rep_all[[col for col in rep_all.columns if 'DNA' in col]].sum(axis=1)

# get number of reads reads
DNA_reads = [rep1DNA, rep2DNA, rep3DNA, rep1DNA+rep2DNA+rep3DNA]
RNA_reads = [rep1RNA, rep2RNA, rep3RNA, rep1RNA+rep2RNA+rep3RNA]
print(DNA_reads)
print(RNA_reads)

# group oligos and sum RNA, DNA counts. Also counts how many barcodes we have per oligo. Calculate ratio and perform log2 on ratio
rep_list = []

for i, rep in enumerate([rep1, rep2, rep3, rep_all]):
    rep = rep[(rep['DNA'] > 0) | (rep['RNA'] > 0)]
    rep[['oligo', 'bc']] = rep.oligo_bc.str.rsplit("_", n=1, expand=True) #06.10.2024 NM - changed  rsplit("_", 1, expand=True) to rsplit("_", n=1, expand=True)
    grouped_df = rep.groupby('oligo')
    final = grouped_df[['RNA', 'DNA']].agg('sum')
    final["count"] = grouped_df.size()
    final["RNA_pseudo"] = final["RNA"] + 1
    final["DNA_pseudo"] = final["DNA"] + 1
    final["RNA_cpm"] = (final["RNA_pseudo"]*1000000)/RNA_reads[i]
    final["DNA_cpm"] = (final["DNA_pseudo"]*1000000)/DNA_reads[i]
    final["ratio"] = final["RNA_cpm"]/final["DNA_cpm"]
    final["ratio_log"] = np.log2(final["ratio"])
    final['oligo_column'] = final.index
    rep_list.append(final[["oligo_column", "ratio_log", 'DNA','RNA',"RNA_pseudo","DNA_pseudo","RNA_cpm","DNA_cpm","ratio","count"]])

rep1_counts = rep_list[0]
rep2_counts = rep_list[1]
rep3_counts = rep_list[2]
rep_all_counts = rep_list[3]

# merge the different outputs
ratio_df = rep1_counts.merge(rep2_counts,on='oligo_column',how="outer",suffixes=("_rep1",None)).merge(rep3_counts,on='oligo_column',how="outer",suffixes=("_rep2",None)).merge(rep_all_counts,on="oligo_column",how="outer",suffixes=("_rep3","_rep_comb"))

ratio_df.to_csv(f'./{cells}/{library}{adaptor}/output/activity_after_filter/ratio_after_filter.csv', header=True, index=False)
