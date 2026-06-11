import pandas as pd
import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]
project = sys.argv[4]

# create list of active sequences (wo archaic/modern)
# do unique on list
# filter input for mpranalyze comparative using this list on locus column


# use quant output instead of direct mpra output
output_quantitative = pd.read_csv(f'./{cells}/{library}{adaptor}/output/activity_after_filter/comb_df_adjusted_fdr.csv')
#output_quantitative = pd.read_csv(f'./{cells}/{library}{adaptor}/output/activity_after_filter/comb_df_adjusted_fdr_L3a3_Hh.csv') # NM 01-11-2024; Toggle for Hh analysis in L3a3

print("Oligo counts before filtering:", len(output_quantitative))
#output_quantitative['oligo'] = output_quantitative['oligo'].str.replace("_ancestral|_derived", '') # remove allele
output_quantitative['oligo'] = output_quantitative['oligo'].str.replace("_derived", '') # remove allele
output_quantitative['oligo'] = output_quantitative['oligo'].str.replace("_ancestral", '') # remove allele


active_df = output_quantitative[output_quantitative["input_comparative"]=="yes"] 
print("Oligo counts after filtering:", len(active_df))

active_oligos = list(set(active_df["oligo"].tolist())) #turn into list using only unique oligos
print("Unique oligo counts:", len(active_oligos))

print(active_df.head().to_string())
print(active_df.info())

# read in input for mpranalyze comparative
input_comparative_rna = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_comparative/{project}/{cells}_{library}{adaptor}_RNA_MPRAnalyze_comparative_improved.txt', sep='\t', header =0)
input_comparative_dna = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_comparative/{project}/{cells}_{library}{adaptor}_DNA_MPRAnalyze_comparative_improved.txt', sep='\t', header =0)

input_comparative_rna_active = input_comparative_rna[input_comparative_rna['locus'].isin(active_oligos)]
input_comparative_dna_active = input_comparative_dna[input_comparative_dna['locus'].isin(active_oligos)]

print(input_comparative_rna_active.info())
print(input_comparative_dna_active.info())

input_comparative_rna_active.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/{project}/{cells}_{library}{adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header=True, index=False)
input_comparative_dna_active.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/{project}/{cells}_{library}{adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr.txt', sep='\t', header=True, index=False)
