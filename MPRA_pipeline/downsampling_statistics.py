#NM 21-11-2024 This code is based on Katharina's downsampling_corr.py.
#The differences: 
#(1) It taks several downsampled libraries at once 
#(2) It produces graphs and tables instead of printing the values 
#(3) It also calculates additional statistics like the number of barcodes. 

import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from scipy import stats
from itertools import islice

library = sys.argv[1]
adaptor = sys.argv[2]

RNA_barcode_list = []
DNA_barcode_list = []
oligo_barcode_list = []

corr_1_2_list = []
corr_2_3_list = []
corr_3_1_list = []

pos_corr_1_2_list = []
pos_corr_2_3_list = []
pos_corr_3_1_list = []

neg_corr_1_2_list = []
neg_corr_2_3_list = []
neg_corr_3_1_list = []

#sactive_oligo_variance_list = []



downsampling_numbers = [1,5,50,70,80,90]
#downsampling_numbers = [1,5]


downsampling_cells = [f"chondrocytesDownsampling{num}" for num in downsampling_numbers]

for cell_type in downsampling_cells:
    print(cell_type)

    UMI_all_reps = pd.read_csv(f'./{cell_type}/{library}{adaptor}/output/UMI/UMI_exploded_std2_filter.txt', sep='\t')
    UMI_all_reps.set_index('oligo_bc', inplace=False)
    rep1DNA=UMI_all_reps[f"DNA_{cell_type}_rep1"].sum()
    rep1RNA=UMI_all_reps[f"RNA_{cell_type}_rep1"].sum()
    rep2DNA=UMI_all_reps[f"DNA_{cell_type}_rep2"].sum()
    rep2RNA=UMI_all_reps[f"RNA_{cell_type}_rep2"].sum()
    rep3DNA=UMI_all_reps[f"DNA_{cell_type}_rep3"].sum()
    rep3RNA=UMI_all_reps[f"RNA_{cell_type}_rep3"].sum()
    
    # print(UMI_all_reps.head())
    # print("###counts########")
    # print("rep1DNA:" ,rep1DNA)
    # print("rep1RNA:" ,rep1RNA)
    # print("rep2DNA:" ,rep2DNA)
    # print("rep2RNA:" ,rep2RNA)
    # print("rep3DNA:" ,rep3DNA)
    # print("rep3RNA:" ,rep3RNA)

    DNA_reads = [rep1DNA, rep2DNA, rep3DNA]
    RNA_reads = [rep1RNA, rep2RNA, rep3RNA]

    rep1 = UMI_all_reps[["oligo_bc", f"RNA_{cell_type}_rep1", f"DNA_{cell_type}_rep1"]]
    rep2 = UMI_all_reps[["oligo_bc", f"RNA_{cell_type}_rep2", f"DNA_{cell_type}_rep2"]]
    rep3 = UMI_all_reps[["oligo_bc", f"RNA_{cell_type}_rep3", f"DNA_{cell_type}_rep3"]]

    rep1 = rep1.rename(columns={f"RNA_{cell_type}_rep1": 'RNA', f"DNA_{cell_type}_rep1": 'DNA'})
    rep2 = rep2.rename(columns={f"RNA_{cell_type}_rep2": 'RNA', f"DNA_{cell_type}_rep2": 'DNA'})
    rep3 = rep3.rename(columns={f"RNA_{cell_type}_rep3": 'RNA', f"DNA_{cell_type}_rep3": 'DNA'})
        
    # Number of unique barcodes in all reps combined
    
    combined_RNA_oligos=[]
    combined_DNA_oligos=[]
    
    for i,rep in enumerate([rep1, rep2, rep3]):
        rep_RNA = rep[(rep["RNA"] > 0)]
        rep_DNA = rep[(rep["DNA"] > 0)]
        combined_RNA_oligos.append(rep_RNA["oligo_bc"])
        combined_DNA_oligos.append(rep_DNA["oligo_bc"])
    

    combined_DNA_oligos = set(pd.concat(combined_DNA_oligos))
    combined_RNA_oligos = set(pd.concat(combined_RNA_oligos))
    
    RNA_barcode_list.append(len(combined_RNA_oligos))
    DNA_barcode_list.append(len(combined_DNA_oligos))
    
    # print("combined_DNA_oligos top 5: ",list(islice(combined_DNA_oligos, 5)))
    # print("combined_RNA_oligos top 5: ",list(islice(combined_RNA_oligos, 5)))

    # Number of unique oligos in all reps combined
    UMI_all_reps_copy = UMI_all_reps.copy()
    UMI_all_reps_copy[['oligo','bc']] = UMI_all_reps_copy.oligo_bc.str.rsplit("_", n=1, expand=True)
    UMI_all_reps_copy = UMI_all_reps_copy[(UMI_all_reps_copy[f"RNA_{cell_type}_rep1"] > 0) | (UMI_all_reps_copy[f"RNA_{cell_type}_rep2"] > 0)| (UMI_all_reps_copy[f"RNA_{cell_type}_rep3"] > 0)| (UMI_all_reps_copy[f"DNA_{cell_type}_rep1"] > 0)| (UMI_all_reps_copy[f"DNA_{cell_type}_rep2"] > 0)| (UMI_all_reps_copy[f"DNA_{cell_type}_rep3"] > 0)]
    all_oligos = set(UMI_all_reps_copy["oligo"])
    oligo_barcode_list.append(len(all_oligos))
    

    # take out barcodes that have RNA and DNA 0
    rep1 = rep1[(rep1["DNA"] > 0) | (rep1["RNA"] > 0)]
    rep2 = rep2[(rep2["DNA"] > 0) | (rep2["RNA"] > 0)]
    rep3 = rep3[(rep3["DNA"] > 0) | (rep3["RNA"] > 0)]

    # group oligos and sum RNA, DNA counts. Also counts how many barcodes we have per oligo. Caluclate ratio and perform log2 on ratio
    rep_list = []

    for i, rep in enumerate([rep1, rep2, rep3]):
        rep[['oligo','bc']] = rep.oligo_bc.str.rsplit("_", n=1, expand=True)
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
        rep_list.append(final)


    rep1_counts = rep_list[0]
    rep2_counts = rep_list[1]
    rep3_counts = rep_list[2]

    # merge different replicates
    rep1_2 = rep1_counts.merge(rep2_counts, how="inner", on="oligo_column")
    rep2_3 = rep2_counts.merge(rep3_counts, how="inner", on="oligo_column")
    rep3_1 = rep3_counts.merge(rep1_counts, how="inner", on="oligo_column")

    # add labels for controls
    rep1_2.loc[rep1_2['oligo_column'].str.contains('NegCtrl_not_active'), 'Label'] = 'Negative control'
    rep1_2.loc[rep1_2['oligo_column'].str.contains(f'PosCtrl_diff_Osteoblasts'), 'Label'] = 'Positive control'
    rep2_3.loc[rep2_3['oligo_column'].str.contains('NegCtrl_not_active'), 'Label'] = 'Negative control'
    rep2_3.loc[rep2_3['oligo_column'].str.contains(f'PosCtrl_diff_Osteoblasts'), 'Label'] = 'Positive control'
    rep3_1.loc[rep3_1['oligo_column'].str.contains('NegCtrl_not_active'), 'Label'] = 'Negative control'
    rep3_1.loc[rep3_1['oligo_column'].str.contains(f'PosCtrl_diff_Osteoblasts'), 'Label'] = 'Positive control'

    rep1_2['Label'] = rep1_2['Label'].fillna(value='No control')
    rep2_3['Label'] = rep2_3['Label'].fillna(value='No control')
    rep3_1['Label'] = rep3_1['Label'].fillna(value='No control')


    rep1_2['Label'] = pd.Categorical(rep1_2.Label, categories=['No control', 'Positive control', 'Negative control'], ordered=True)
    rep1_2.sort_values('Label', inplace = True)
    rep2_3['Label'] = pd.Categorical(rep2_3.Label, categories=['No control', 'Positive control', 'Negative control'], ordered=True)
    rep2_3.sort_values('Label', inplace = True)
    rep3_1['Label'] = pd.Categorical(rep3_1.Label, categories=['No control', 'Positive control', 'Negative control'], ordered=True)
    rep3_1.sort_values('Label', inplace = True)

    # filter to take only oligos with specific amount of barcodes
    rep10bcs = []
    rep5bcs = []
    rep_all_bcs = []

    # print(f"Length of rep1_2: {len(rep1_2)}")

    # calculate overall correlation between replicates >=10 bcs or >=5 bcs

    corr_1_2_log = rep1_2["ratio_log_x"].corr(rep1_2["ratio_log_y"])
    corr_2_3_log = rep2_3["ratio_log_x"].corr(rep2_3["ratio_log_y"])
    corr_3_1_log = rep3_1["ratio_log_x"].corr(rep3_1["ratio_log_y"])
    
    corr_1_2_list.append(corr_1_2_log)
    corr_2_3_list.append(corr_2_3_log)
    corr_3_1_list.append(corr_3_1_log)
    
    # calculate positive controls
    
    pos_1_2 = rep1_2[rep1_2["Label"] == "Positive control"]
    pos_2_3 = rep2_3[rep2_3["Label"] == "Positive control"]
    pos_3_1 = rep3_1[rep3_1["Label"] == "Positive control"]
    
    pos_corr_1_2_log = pos_1_2["ratio_log_x"].corr(pos_1_2["ratio_log_y"])
    pos_corr_2_3_log = pos_2_3["ratio_log_x"].corr(pos_2_3["ratio_log_y"])
    pos_corr_3_1_log = pos_3_1["ratio_log_x"].corr(pos_3_1["ratio_log_y"])    

    pos_corr_1_2_list.append(pos_corr_1_2_log)
    pos_corr_2_3_list.append(pos_corr_2_3_log)
    pos_corr_3_1_list.append(pos_corr_3_1_log)
    
    if cell_type == "chondrocytesDownsampling50":
        datasets = [
            ("1-2", pos_1_2["ratio_log_x"], pos_1_2["ratio_log_y"]),
            ("2-3", pos_2_3["ratio_log_x"], pos_2_3["ratio_log_y"]),
            ("3-1", pos_3_1["ratio_log_x"], pos_3_1["ratio_log_y"]),
        ]
    
        plt.clf()
    
        fig, axes = plt.subplots(3, 1, figsize=(8, 12))
    
        for i, (label, ratio_log_x, ratio_log_y) in enumerate(datasets):
            ax = axes[i]
            ax.scatter(ratio_log_x, ratio_log_y, alpha=0.7)
            ax.set_title(f"{label} Scatter Plot")
            ax.set_xlabel("Ratio Log X")
            ax.set_ylabel("Ratio Log Y")
    
        # Adjust layout
        plt.tight_layout()
    
        # Save the figure
        plt.savefig(f"./additional/downsampling_corr_barcodes_oligos_Nadav/L4a1/correlation_scatterplots_{cell_type}.png", dpi=1000)
    
    # calculate negative controls
    neg_1_2 = rep1_2[rep1_2["Label"] == "Negative control"]
    neg_2_3 = rep2_3[rep2_3["Label"] == "Negative control"]
    neg_3_1 = rep3_1[rep3_1["Label"] == "Negative control"]
    
    neg_corr_1_2_log = neg_1_2["ratio_log_x"].corr(neg_1_2["ratio_log_y"])
    neg_corr_2_3_log = neg_2_3["ratio_log_x"].corr(neg_2_3["ratio_log_y"])
    neg_corr_3_1_log = neg_3_1["ratio_log_x"].corr(neg_3_1["ratio_log_y"])    
    
    neg_corr_1_2_list.append(neg_corr_1_2_log)
    neg_corr_2_3_list.append(neg_corr_2_3_log)
    neg_corr_3_1_list.append(neg_corr_3_1_log) 


    # calculate variance of active oligos
    
    
    active_oligo_variance_list = []

    
        

print("Done!")

#Plot the correlations

plt.clf()
#All
sns.lineplot(x=downsampling_numbers,y=corr_1_2_list, marker="o",color="blue", linestyle="-",label = "1-2")
sns.lineplot(x=downsampling_numbers,y=corr_2_3_list, marker="o",color="yellow", linestyle="-",label = "2-3")
sns.lineplot(x=downsampling_numbers,y=corr_3_1_list, marker="o",color="green", linestyle="-",label = "1-3")
#Positive
sns.lineplot(x=downsampling_numbers,y=pos_corr_1_2_list, marker="P",color="blue", linestyle="-",label = "1-2 Positive controls")
sns.lineplot(x=downsampling_numbers,y=pos_corr_2_3_list, marker="P",color="yellow", linestyle="-",label = "2-3 Positive controls")
sns.lineplot(x=downsampling_numbers,y=pos_corr_3_1_list, marker="P",color="green", linestyle="-",label = "1-3 Positive controls")
#Negative
sns.lineplot(x=downsampling_numbers,y=neg_corr_1_2_list, marker="^",color="blue", linestyle="-",label = "1-2 Negative controls")
sns.lineplot(x=downsampling_numbers,y=neg_corr_2_3_list, marker="^",color="yellow", linestyle="-",label = "2-3 Negative controls")
sns.lineplot(x=downsampling_numbers,y=neg_corr_3_1_list, marker="^",color="green", linestyle="-",label = "1-3 Negative controls")
plt.xlabel("% downsampling")
plt.ylabel("Correlation")
plt.xticks(ticks=downsampling_numbers)
plt.title("Barcode counts")
plt.ylim(0,1)

plt.savefig("./additional/downsampling_corr_barcodes_oligos_Nadav/L4a1/correlations.png", dpi=1000)


#Plot the barcode counts
plt.clf()
sns.lineplot(x=downsampling_numbers,y=RNA_barcode_list,color="blue", linestyle="-",label = "RNA")
sns.lineplot(x=downsampling_numbers,y=DNA_barcode_list,color="yellow", linestyle="-",label = "DNA")
plt.xlabel("% downsampling")
plt.ylabel("Barcode counts")
plt.xticks(ticks=downsampling_numbers)
plt.title("Oligos counts")
plt.savefig("./additional/downsampling_corr_barcodes_oligos_Nadav/L4a1/barcode_counts.png", dpi=1000)

#Plot the oligo counts
plt.clf()
sns.lineplot(x=downsampling_numbers,y=oligo_barcode_list, linestyle="-")
plt.xlabel("% downsampling")
plt.ylabel("Oligo counts")
plt.xticks(ticks=downsampling_numbers)
plt.title("Correlations between replicated")
plt.tight_layout()
plt.savefig("./additional/downsampling_corr_barcodes_oligos_Nadav/L4a1/oligo_counts.png", dpi=1000)


