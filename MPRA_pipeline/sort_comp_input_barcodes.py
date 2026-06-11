import pandas as pd
import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

# set the follwing to TRUE if running only on additional sequences that were added after running on fdr on combined libraries:

if False:
    additional_sequences = "_additional_sequences"
else:
        additional_sequences = ""

input_comparative_dna = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr{additional_sequences}.txt', sep='\t', header =0)
input_comparative_rna = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr{additional_sequences}.txt', sep='\t', header =0)

# counting bcs by countung non-zeros
input_comparative_rna["count_bcs_rna"] = (input_comparative_rna != 0).sum(axis=1)
input_comparative_dna["count_bcs_dna"] = (input_comparative_dna != 0).sum(axis=1)

input_comparative_dna["count_bcs_rna"] = input_comparative_rna["count_bcs_rna"]

# sort both dataframes according to rna bc counts
input_comparative_rna.sort_values(by='count_bcs_rna', ascending=False, inplace = True)
input_comparative_dna.sort_values(by='count_bcs_rna', ascending=False, inplace = True)


print(input_comparative_rna[['locus','count_bcs_rna']].head(30))
print(input_comparative_dna[['locus','count_bcs_dna', 'count_bcs_rna']].head(10))

input_comparative_rna_to_save = input_comparative_rna.drop(columns=['count_bcs_rna'])
input_comparative_dna_to_save = input_comparative_dna.drop(columns=['count_bcs_dna', 'count_bcs_rna'])

# save sorted wihtout column headers - later used as input for csplit
input_comparative_rna_to_save.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted{additional_sequences}.txt', sep='\t', header=False, index=False)
input_comparative_dna_to_save.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/{cells}_{library}{adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted{additional_sequences}.txt', sep='\t', header=False, index=False)

# loop over count bcs rna to get split positions

count_bcs_rna_list = input_comparative_rna['count_bcs_rna'].tolist()
chunk_sum = 0
split_indices = []
# max_bc= max(count_bcs_rna_list) # does not work well because some libraries have high max bc count - chunks will be to large
max_bc= 8000 # 8000 seems to be reasonable (at least in humanMPRA)
print(max_bc)

for idx, count_bcs in enumerate(count_bcs_rna_list):
    if chunk_sum + count_bcs > max_bc:
        if idx != 0: # first cannot be 0 because then the first chunk will be empty
            split_indices.append(idx)
        chunk_sum = 0
    chunk_sum += count_bcs

split_indices_plus_one = [x+1 for x in split_indices] # have to add 1 because of how csplit uses the positions
print(len(split_indices_plus_one))


with open(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_comparative/split_positions{additional_sequences}.txt', 'w') as f:
    for line in split_indices_plus_one:
        f.write(f"{line} ")