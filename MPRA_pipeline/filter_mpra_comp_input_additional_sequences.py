# this is specifically for the chondrocytes. because of achange in how we do the fdr the list of oligos we want to include in comparative changed
# here we take the the oligos that now also have to be included and create another input for mpranalyze comparative

import pandas as pd
import sys     
      
cells = sys.argv[1]

# read in the combined output from mpranalyze quantitative
output_quantitative = pd.read_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv')

output_quantitative['oligo'] = output_quantitative['oligo'].str.replace("_ancestral", '') # remove allele
output_quantitative['oligo'] = output_quantitative['oligo'].str.replace("_derived", '') # remove allele

# create a list of all oligos that were defined as not active before and are now defined as active
new_active_df = output_quantitative[(output_quantitative['activity_adjusted_combined'] == 'active') & (output_quantitative['activity_adjusted'] != 'active')] 
new_active_oligos = list(set(new_active_df["oligo"].tolist())) #turn into list using only unique oligos

# also create a list of oligos that were previously defined as active
old_active_df = output_quantitative[output_quantitative["input_comparative"]=="yes"] 
old_active_oligos = list(set(old_active_df["oligo"].tolist())) #turn into list using only unique oligos

# make sure that the new active oligo list does not contain oligos from the old list - this could happen because the oligos are in pairs
# one part of the pair might be defined as active only now but when the other member was already defined as active previously we already included it before and do not have to include it again
new_active_oligos_updated = [item for item in new_active_oligos if item not in old_active_oligos]
print(len(new_active_oligos))
print(len(old_active_oligos))
print(len(new_active_oligos_updated))

# loop of the libraries and creata new input files for the new subset
for library in ["L1", "L2", "L3","L4"]:
    for adaptor in ["a1", "a2", "a3"]:
        if library=="L4" and adaptor in("a2","a3"):
            continue
        input_comparative_rna = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_RNA_MPRAnalyze_comparative_improved.txt', sep='\t', header =0)
        input_comparative_dna = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_DNA_MPRAnalyze_comparative_improved.txt', sep='\t', header =0)

        input_comparative_rna_active = input_comparative_rna[input_comparative_rna['locus'].isin(new_active_oligos_updated)]
        input_comparative_dna_active = input_comparative_dna[input_comparative_dna['locus'].isin(new_active_oligos_updated)]

        input_comparative_rna_active.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_additional_sequences.txt', sep='\t', header=True, index=False)
        input_comparative_dna_active.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_additional_sequences.txt', sep='\t', header=True, index=False)
