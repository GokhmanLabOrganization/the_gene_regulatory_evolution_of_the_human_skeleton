# chnages to incorporate:
# use normalized rna/dna counts
# use mad score instead of alpha x
# replace na in diff activity fdr with 1
# also in fdr for activity?
# change the way sequences are added - right now, if the oligo did not have association the sequence does not appear x

import pandas as pd
import sys
from Bio import SeqIO

cells = sys.argv[1]

quant = pd.read_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv', header=0, usecols=["oligo","ratio_log_rep_comb", "DNA_rep_comb", "RNA_rep_comb","count_rep_comb","alpha", "mad.score", "fdr.mad_adjusted_combined", "activity_adjusted_combined", "control_type",'orientation_fix'])
print("number of rows in original DF:",quant.shape[0])
quant.rename(columns = {'ratio_log_rep_comb':'RNA_DNA_ratio_log2', 'DNA_rep_comb':'DNA_counts_raw', 'RNA_rep_comb':'RNA_counts_raw', 'count_rep_comb':'barcode_count', 'alpha':'alpha_activity_estimate', 'mad.score':"normalized_activity_estimate", 'fdr.mad_adjusted_combined':'activity_fdr', 'activity_adjusted_combined':'activity'}, inplace = True)
quant = quant.loc[~(quant['orientation_fix']== 'fixed_in_L4')] #NM 18-12-2024 to remove oligos with fixed orientation in L4a1
quant.drop(columns=['orientation_fix'], inplace=True)

print("after removing oligos which fixed orientation in L4:", quant.shape[0])

#remove controls and non-relevant oligos
quant = quant[quant.control_type == "No control"]
print("after removing controls: ",quant.shape[0])

quant = quant[~quant['oligo'].str.contains('_Hh_')]
quant = quant[~quant['oligo'].str.contains('hh.SCREEN')]
quant = quant[~quant['oligo'].str.contains('FABP7')]
quant = quant[~quant['oligo'].str.contains('hh.missing.oligos')]
quant.drop(['control_type'], axis=1, inplace=True)
print("after removing HH and non-relevant oligos: ",quant.shape[0])


# add normalized counts from ratio df


# Add a new column for column suffix (ancestral or derived)
quant['allele'] = quant['oligo'].apply(lambda x: 'ancestral' if 'ancestral' in x else 'derived' if 'derived' in x else None)

#remove allele from oligo
quant['oligo'] = quant['oligo'].str.replace("_ancestral|_derived", '',regex=True)


# Pivot the DataFrame
quant_wide = quant.pivot(index='oligo', columns='allele', values=['RNA_DNA_ratio_log2', 'DNA_counts_raw', 'RNA_counts_raw', 'barcode_count', 'alpha_activity_estimate', 'normalized_activity_estimate', 'activity_fdr', 'activity']).reset_index()
print(quant_wide.shape[0])
# Flatten the multi-level columns
quant_wide.columns = [f'{col[0]}_{col[1]}' if col[1] else col[0] for col in quant_wide.columns]
print(quant_wide.shape[0])
# print head
# print(quant_wide.head().to_string())

# merge with comp output
comp = pd.read_csv(f'./{cells}/comparative_analysis_combined/mpranalyze_comp_res_fdr_complete_w_add_seq.csv', header=0, usecols=["oligo","logFC","fdr.wo_controls", "differntial.wo_controls"])
comp.rename(columns = {'logFC':'logFC_derived_vs_ancestral', 'fdr.wo_controls':'differential_activity_fdr', 'differntial.wo_controls':'differential_activity'}, inplace = True)

quant_comp = quant_wide.merge(comp, how="left", on="oligo")
print(quant_comp.shape[0])
# remove Hh
quant_comp = quant_comp[~quant_comp['oligo'].str.contains('_Hh_')]
quant_comp = quant_comp[~quant_comp['oligo'].str.contains('hh.SCREEN')]
quant_comp = quant_comp[~quant_comp['oligo'].str.contains('FABP7')]
quant_comp = quant_comp[~quant_comp['oligo'].str.contains('hh.missing.oligos')]

print("Counts after removing HH only oligos", quant_comp.shape[0])

# add sequences
# Initialize empty lists to store modified identifiers and sequences
identifiers = []
sequences = []

# Parse the FASTA file
fasta_file = f"./oligo_fasta/L1L2L3L4a1.fasta" # NM 16-12-2024 Added L4a1
for record in SeqIO.parse(fasta_file, "fasta"):
    identifier = record.id  # Remove the first character
    identifiers.append(identifier)
    sequence = str(record.seq)
    sequences.append(sequence)

sequence_df = pd.DataFrame({'oligo': identifiers, 'sequence': sequences})
sequence_df['sequence'] = sequence_df['sequence'].str[15:-15]

# Add a new column for column suffix (ancestral or derived)
sequence_df['allele'] = sequence_df['oligo'].apply(lambda x: 'ancestral' if 'ancestral' in x else 'derived' if 'derived' in x else None)

#remove allele from oligo
sequence_df['oligo'] = sequence_df['oligo'].str.replace("_ancestral|_derived", '', regex=True)

#remove some controls which were not removed so far:
sequence_df.dropna(subset=['allele'], inplace=True) #Not sure why I still need to do this

# Pivot the DataFrame
print(sequence_df.loc[sequence_df["allele"].isna()])

sequence_df_wide = sequence_df.pivot(index='oligo', columns='allele', values=["sequence"]).reset_index()
print(sequence_df_wide.shape[0])
# Flatten the multi-level columns
sequence_df_wide.columns = [f'{col[0]}_{col[1]}' if col[1] else col[0] for col in sequence_df_wide.columns]
print(sequence_df_wide.shape[0])

sequence_df_wide = sequence_df_wide[["oligo", "sequence_ancestral", "sequence_derived"]]
# merge sequences to quant_Comp
quant_comp = quant_comp.merge(sequence_df_wide, on="oligo", how="left")

# save df
quant_comp.to_csv(f'./{cells}/comparative_analysis_combined/humanMPRA_with_seq_final2.csv', header=True, index = False)
