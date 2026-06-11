import pandas as pd
import sys
import statsmodels.stats.multitest as smm
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import seaborn as sns
import os

cells = sys.argv[1]

data_df_path = f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv' # NM 14-11-2024

def create_comb_df_combined(cells):
    # read in comb df (fdr adjusted) 
    comb_df_list = []
    for library_adaptor in ["L1a1","L1a2","L1a3","L2a1","L2a2","L2a3","L3a1","L3a2","L3a3","L4a1"]:
        comb_df_list.append(pd.read_csv(f'./{cells}/{library_adaptor}/output/activity_after_filter/comb_df_adjusted_fdr.csv', header=0))

    # concatenate comb dfs
    comb_df_concat = pd.concat(comb_df_list)
    
    print("# oligos before filtering:",len(comb_df_concat))
    #filtered_comb_df_concat = comb_df_concat[(comb_df_concat['DNA_rep_comb'] >= 5) & (comb_df_concat['control_type'] == 'No control')& (~comb_df_concat['pval.mad'].isna())& (~comb_df_concat['oligo'].str.contains('_Hh_'))]
    
    # add col of genomic coords
    comb_df_concat["coords"] = comb_df_concat["oligo"].str.extract(r'(chr[\w\d]+:\d+-\d+)')
    
        
    oligos_with_fixed_orientation = comb_df_concat.loc[comb_df_concat["oligo"].str.contains("promoter.orientation.fix"),"coords"]
    print("number of oligos with their orientation fixed in L4a1:",len(oligos_with_fixed_orientation))
    
    print(oligos_with_fixed_orientation[0:10])
    
    
    #add a new col for test oligos which had their orientation fixed in L4, to be removed from the FDR 
    comb_df_concat["orientation_fix"] = "unchanged"
    comb_df_concat.loc[(
        (comb_df_concat["coords"].isin(oligos_with_fixed_orientation)) &
        (~comb_df_concat["oligo"].str.contains("a1_L4")) &
        (~comb_df_concat["oligo"].str.contains("Ctrl")) &
        (~comb_df_concat["oligo"].str.contains("scrambled"))),
        "orientation_fix"
    ] = "fixed_in_L4"
    
    mask = (
        (comb_df_concat['DNA_rep_comb'] >= 5) & 
        (comb_df_concat['control_type'] == 'No control')& 
        (~comb_df_concat['pval.mad'].isna())& 
        (~comb_df_concat['oligo'].str.contains('_Hh_'))& # Relevant to L3a3 and L4a1
        (~comb_df_concat['oligo'].str.contains('FABP7'))& # Relevant to L4a1
        (~comb_df_concat['oligo'].str.contains('hh.SCREEN'))& # Relevant to L4a1
        (~comb_df_concat['oligo'].str.contains('hh.missing.oligos'))& # Relevant to L4a1
        (~(comb_df_concat['orientation_fix'] == "fixed_in_L4"))
    )
    
    # Add a col named - orientaion fixed in L4 (yes/no) - yes - not in L4 and has locs which are in L4
    # add to mask - remove all oligos which are fixed
    
    filtered_comb_df_concat = comb_df_concat[mask]
    print("# oligos after filtering:",len(filtered_comb_df_concat))
    print(len(filtered_comb_df_concat.index))
    
    rej, pval_corr = smm.fdrcorrection(filtered_comb_df_concat["pval.mad"])
    print(len(pval_corr))
    comb_df_concat.loc[mask , 'fdr.mad_adjusted_combined'] = pval_corr
    comb_df_concat['fdr.mad_adjusted_combined'].fillna(1, inplace=True)

    # add new activity columns
    comb_df_concat.loc[comb_df_concat["fdr.mad_adjusted_combined"] < 0.05, "activity_adjusted_combined"] = "active"
    activity_mask = mask &(comb_df_concat["fdr.mad_adjusted_combined"] >= 0.05)
    comb_df_concat.loc[activity_mask,'activity_adjusted_combined'] = "non_active"
    
    # comb_df_concat.loc[(comb_df_concat['DNA_rep_comb'] >= 5) & 
        # (comb_df_concat['control_type'] == 'No control')& 
        # (~comb_df_concat['pval.mad'].isna()) & 
        # (comb_df_concat["fdr.mad_adjusted_combined"] >= 0.05)& 
        # (~comb_df_concat['oligo'].str.contains('_Hh_')), 'activity_adjusted_combined'] = "non_active"

    # add column whether oligo is input or not for mpranalyze comp
    comb_df_concat["input_comparative_combined"] = "no"
    comb_df_concat.loc[comb_df_concat["fdr.mad_adjusted_combined"] < 0.05, "input_comparative_combined"] = "yes"
    comb_df_concat.loc[comb_df_concat['control_type'] == "PosCtrl_diff", "input_comparative_combined"] = "yes"

    # count how many were defined as active when looking at libraries separately and are not defined active when looking at libraries combined (and the other way around)
    print("active in combined, not active in separately")
    print(len(comb_df_concat[(comb_df_concat['activity_adjusted_combined'] == 'active') & (comb_df_concat['activity_adjusted'] != 'active')]))
    print("active separately, not active in combined")
    print(len(comb_df_concat[(comb_df_concat['activity_adjusted'] == 'active') & (comb_df_concat['activity_adjusted_combined'] != 'active')]))

    # count how many were input for comparative when looking at libraries separately and are not input for comparative when looking at libraries combined (and the other way around)
    print("input in combined, not input separately")
    print(len(comb_df_concat[(comb_df_concat['input_comparative_combined'] == 'yes') & (comb_df_concat['input_comparative'] != 'yes')]))
    print("input separately, not input in combined")
    print(len(comb_df_concat[(comb_df_concat['input_comparative'] == 'yes') & (comb_df_concat['input_comparative_combined'] != 'yes')]))

    #save df
    comb_df_concat.to_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv', header=True, index = False)
    small_df = comb_df_concat[['oligo', 'coords','mad.score', 'pval.mad', 'count_rep_comb', 'ratio_log_rep_comb', 'fdr.mad_adjusted', 'fdr.mad_adjusted_combined', 'activity_adjusted', 'activity_adjusted_combined', 'input_comparative', 'input_comparative_combined',"orientation_fix"]]
    small_df.to_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr_small.csv', header=True, index = False)
    
    return comb_df_concat


def main():
    if not os.path.exists(data_df_path):
        comb_df_concat = create_comb_df_combined(cells)
        return comb_df_concat
        print('comb_df_combined was created')
    else:
        comb_df_concat = pd.read_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv')
        return comb_df_concat
        print('used exisiting comb_df_combined')

comb_df_concat = main()

num_bins = 20

#check whether mad score normalizes well between libraries by looking at the distribution of shared sequences
positive_control_lib_subset = []
negative_control_lib_subset = []
for library_adaptor in ["a1_L1","a2_L1","a3_L1","a1_L2","a2_L2","a3_L2","a1_L3","a2_L3","a3_L3","a1_L4"]:
    mad_positive = comb_df_concat.loc[comb_df_concat['oligo'].str.contains(f'{library_adaptor}') & comb_df_concat['oligo'].str.contains('PosCtrl'), 'mad.score']
    mad_negative = comb_df_concat.loc[comb_df_concat['oligo'].str.contains(f'{library_adaptor}') & comb_df_concat['oligo'].str.contains('NegCtrl_not_active'), 'mad.score']
    positive_control_lib_subset.append(mad_positive)
    negative_control_lib_subset.append(mad_negative) 
plt.clf()
plt.hist(positive_control_lib_subset, label=["L1a1","L1a2","L1a3","L2a1","L2a2","L2a3","L3a1","L3a2","L3a3","L4a1"], bins=20, histtype='step', fill=False)
plt.legend(loc='upper right')
plt.savefig(f'./{cells}/quantitative_analysis_combined/hist_mad_score_positive_controls.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/hist_mad_score_positive_controls.png', dpi=330)
plt.clf()
plt.hist(negative_control_lib_subset, label=["L1a1","L1a2","L1a3","L2a1","L2a2","L2a3","L3a1","L3a2","L3a3","L4a1"], histtype='step', fill=False)
plt.legend(loc='upper right')
plt.savefig(f'./{cells}/quantitative_analysis_combined/hist_mad_score_negative_controls.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/hist_mad_score_negative_controls.png', dpi=330)


plt.clf()
fig, axes = plt.subplots(10, 2, sharex='col', figsize=(7, 14))
for n, library_adaptor in enumerate(["a1_L1","a2_L1","a3_L1","a1_L2","a2_L2","a3_L2","a1_L3","a2_L3","a3_L3","a1_L4"]):
    mad_positive = comb_df_concat.loc[comb_df_concat['oligo'].str.contains(f'{library_adaptor}') & (comb_df_concat['oligo'].str.contains('PosCtrl_osteo')|comb_df_concat['oligo'].str.contains('PosCtrl_neuron')|comb_df_concat['oligo'].str.contains('PosCtrl_chondro'))]
    mad_negative = comb_df_concat.loc[comb_df_concat['oligo'].str.contains(f'{library_adaptor}') & comb_df_concat['oligo'].str.contains('NegCtrl_not_active')]
    sns.histplot(ax=axes[n, 0], data=mad_positive, x=f"mad.score",bins = num_bins)
    sns.histplot(ax=axes[n, 1], data=mad_negative, x=f"mad.score",bins = num_bins)
    axes[n, 0].set_title(f'{library_adaptor}, pos')
    axes[n, 1].set_title(f'{library_adaptor}, neg')
plt.suptitle("madscore distribution for positive and negative controls")
plt.savefig(f'./{cells}/quantitative_analysis_combined/hist_mad_score.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/hist_mad_score.png', dpi=330)


plt.clf()
sns.histplot(data = comb_df_concat,x= "fdr.mad_adjusted_combined",bins=20)
plt.title("Histogram of FDR adjusted P value for all libraries cobined")
plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_hist_p_val.png', dpi=330)



plt.clf()
fig, axes = plt.subplots(5, 2, sharex='col', figsize=(10, 20), sharey=True)

for n, library_adaptor in enumerate(["a1_L1","a2_L1","a3_L1","a1_L2","a2_L2","a3_L2","a1_L3","a2_L3","a3_L3","a1_L4"]):
    row, col = divmod(n, 2)
    sns.histplot(ax=axes[row, col], data=comb_df_concat.loc[comb_df_concat['oligo'].str.contains(f'{library_adaptor}')], x=f"pval.mad",bins = num_bins)
    axes[row,col].set_title(library_adaptor)

plt.tight_layout()

plt.suptitle("P val Distribution Per Library")
plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_hist_pval_per_lib.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_hist_pval_per_lib.png', dpi=330)

#Pie charts
plt.clf()

libs = ["a1_L1", "a2_L1", "a3_L1", "a1_L2", "a2_L2", "a3_L2", "a1_L3", "a2_L3", "a3_L3", "a1_L4"]
counts = {lib: comb_df_concat['oligo'].str.contains(lib).sum() for lib in libs}

plt.figure(figsize=(8, 8))
plt.pie(counts.values(), labels=counts.keys(), autopct='%1.1f%%', startangle=90)
plt.title('All oligos divided by library')
plt.show()

plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_pie_all_oligos.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_pie_all_oligos.png', dpi=330)

plt.clf()

input_comp_df = comb_df_concat[comb_df_concat["input_comparative"] == "yes"]

counts = {lib: input_comp_df['oligo'].str.contains(lib).sum() for lib in libs}

plt.figure(figsize=(8, 8))
plt.pie(counts.values(), labels=counts.keys(), autopct='%1.1f%%', startangle=90)
plt.title('input_comparative oligos divided by library')
plt.show()

plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_pie_input_comparative_oligos.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_pie_input_comparative_oligos.png', dpi=330)

plt.clf()

input_comp_df = comb_df_concat[comb_df_concat["input_comparative_combined"] == "yes"]
counts = {lib: input_comp_df['oligo'].str.contains(lib).sum() for lib in libs}

plt.figure(figsize=(8, 8))
plt.pie(counts.values(), labels=counts.keys(), autopct='%1.1f%%', startangle=90)
plt.title('input_comparative_combined oligos divided by library')
plt.show()

plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_pie_input_comparative_combined_oligos.pdf')
plt.savefig(f'./{cells}/quantitative_analysis_combined/{cells}_pie_input_comparative_combined_oligos.png', dpi=330)