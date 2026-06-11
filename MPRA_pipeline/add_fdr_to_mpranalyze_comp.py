import pandas as pd
import sys
import statsmodels.stats.multitest as smm
import os

cells = sys.argv[1]

include_additional_sequences = True
# read in concataneted output of mpranalyze comparative for all libraries
output_list = []
# for library_adaptor in ["L1a1","L1a2","L1a3","L2a1","L2a2","L2a3","L3a1","L3a2","L3a3","L4a1"]:
    # path = f'./{cells}/{library_adaptor}/output/mpranalyze_comparative/'
    # output_list.append(pd.read_csv(path+'mpranalyze_comp_res_filter_sorted.txt', sep='\t', header=0, names = ['oligo','statistic','pval','fdr','df.test','df.dna','df.rna.full','df.rna.red','logFC']))
    # print("added original oligos of",library_adaptor)
    # if include_additional_sequences:
        # if os.path.exists(path+'additional_sequences_mpranalyze_comp_res_filter_sorted.txt'):
            # output_list.append(pd.read_csv(path+'additional_sequences_mpranalyze_comp_res_filter_sorted.txt', sep='\t', header=0, names = ['oligo','statistic','pval','fdr','df.test','df.dna','df.rna.full','df.rna.red','logFC']))
            # print("added additional oligos of",library_adaptor)

for library_adaptor in ["L1a1","L1a2","L1a3","L2a1","L2a2","L2a3","L3a1","L3a2","L3a3","L4a1"]:
    path = f'./{cells}/{library_adaptor}/output/mpranalyze_comparative/'
    original_seq = pd.read_csv(path+'mpranalyze_comp_res_filter_sorted.txt', sep='\t', header=0, names = ['oligo','statistic','pval','fdr','df.test','df.dna','df.rna.full','df.rna.red','logFC'])
    print("added original oligos of",library_adaptor)
    og_olig_count = len(original_seq)
    original_seq = original_seq.drop_duplicates(subset='oligo', keep='first')
    print("added original oligos of",library_adaptor)
    og_olig_wo_dup_count = len(original_seq)
    print("removed", og_olig_count-og_olig_wo_dup_count," duplicated original oligos")

    og_olig_count = len(original_seq)
    if include_additional_sequences:
        if os.path.exists(path+'additional_sequences_mpranalyze_comp_res_filter_sorted.txt'):
            additional_seq = pd.read_csv(path+'additional_sequences_mpranalyze_comp_res_filter_sorted.txt', sep='\t', header=0, names = ['oligo','statistic','pval','fdr','df.test','df.dna','df.rna.full','df.rna.red','logFC'])
            add_olig_count = len(additional_seq)
            print("added additional oligos of",library_adaptor)
            
            additional_seq = additional_seq.drop_duplicates(subset='oligo', keep='first')
            add_olig_count_after_dup_rm = len(additional_seq)
            print("removed", add_olig_count-add_olig_count_after_dup_rm," additional oligos which were duplicates of other additional oligos")

            original_seq = original_seq[~original_seq['oligo'].isin(additional_seq['oligo'])]
            og_olig_count_filtered = len(original_seq)
            print("removed", og_olig_wo_dup_count-og_olig_count_filtered," original oligos which were suppoused to be additional")
            output_list.append(original_seq)
            output_list.append(additional_seq)
        else:
            output_list.append(original_seq)
    

 
    
        
        
print(output_list[0].dtypes)
print(output_list[0].head().to_string())
# concatenate outputs
output = pd.concat(output_list)

#TAKE OUT ROWS WHERE WE DONT HAVE FOLDCHANGE
print(output['logFC'].isna().sum())
output = output.dropna(subset=['logFC'])
print(len(output.index))

# take out rows that were previously defined as active but now not anymore - basically only take rows that are now defined as active
# create set with sequences defined as active according to new calculation
output_quant_combined = pd.read_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv')
output_quant_combined['oligo'] = output_quant_combined['oligo'].str.replace("_ancestral", '')
output_quant_combined['oligo'] = output_quant_combined['oligo'].str.replace("_derived", '')

active_df = output_quant_combined[(output_quant_combined['activity_adjusted_combined'] == 'active')] 
active_oligos = list(set(active_df["oligo"].tolist())) #turn into list using only unique oligos

output = output[output['oligo'].isin(active_oligos)]
print(len(output.index))

# do correction for all
rej, pval_corr = smm.fdrcorrection(output["pval"])
output["differential.full"] = rej
output["fdr.full"] = pval_corr

# do correlation for all expect controls
output_wo_controls = output[~output['oligo'].str.contains('Ctrl')]
print(len(output_wo_controls.index))
rej, pval_corr = smm.fdrcorrection(output_wo_controls["pval"])
print(len(pval_corr))
output.loc[~output['oligo'].str.contains('Ctrl'), 'fdr.wo_controls'] = pval_corr
output.loc[~output['oligo'].str.contains('Ctrl'), 'differntial.wo_controls'] = rej

# count differntial/non differential
num_differential = (output["differntial.wo_controls"] == True).sum()
num_not_differential = (output["differntial.wo_controls"] == False).sum()
print(f'# differential: {num_differential}')
print(f'# not differential: {num_not_differential}')
print(f'# total: {num_differential+num_not_differential}')
print(f'% differential: {num_differential/(num_differential+num_not_differential)}')

# save new df
if include_additional_sequences:
    output.to_csv(f'./{cells}/comparative_analysis_combined/mpranalyze_comp_res_fdr_complete_w_add_seq.csv', header=True, index = False)
else:
    output.to_csv(f'./{cells}/comparative_analysis_combined/mpranalyze_comp_res_fdr_complete.csv', header=True, index = False)
