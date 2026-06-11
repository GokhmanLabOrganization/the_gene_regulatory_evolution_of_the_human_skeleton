# create bad for fasta

import pandas as pd
import sys
import re


# get fasta ids and save as df
data = []
with open(f'./oligo_fasta/L1L2L3L4a1.fasta') as fp: #NM 16-12-2024 added L4a1
    for line in fp:
        if line.startswith('>'):
            ident = line[1:].strip()
            data.append(ident)
df = pd.DataFrame(data, columns=['oligo_column'])

# filter out controls
df.loc[df['oligo_column'].str.contains('PosCtrl_osteoblast_active'), 'control_type'] = 'PosCtrl_osteoblast_active'
df.loc[df['oligo_column'].str.contains('PosCtrl_chondrocyte_active'), 'control_type'] = 'PosCtrl_chondrocyte_active'
df.loc[df['oligo_column'].str.contains('PosCtrl_neuron_active'), 'control_type'] = 'PosCtrl_neuron_active'
df.loc[df['oligo_column'].str.contains('NegCtrl_non_SCREEN'), 'control_type'] = 'NegCtrl_non_SCREEN'
df.loc[df['oligo_column'].str.contains('NegCtrl_active_not_diff'), 'control_type'] = 'NegCtrl_active_not_diff'
df.loc[df['oligo_column'].str.contains('scrambled_control'), 'control_type'] = 'scrambled_control'
df.loc[df['oligo_column'].str.contains('NegCtrl_not_active'), 'control_type'] = 'NegCtrl_not_active'
df.loc[df['oligo_column'].str.contains('PosCtrl_diff'), 'control_type'] = 'PosCtrl_diff'
df['control_type'] = df['control_type'].fillna(value='No control')

df_without_controls = df[df["control_type"] == "No control"]
print(df_without_controls.shape[0])

# remove allele from oligo name
df_without_controls['oligo_column'] = df_without_controls['oligo_column'].str.replace("_derived|_ancestral", '',regex=True)
#df_without_controls['oligo_column'] = df_without_controls['oligo_column'].str.replace("_ancestral", '')

print(df_without_controls.shape[0])

# remove duplicates based on oligo column (ancestral/derived)
df_without_controls.drop_duplicates(subset=['oligo_column'], keep = 'first', inplace=True)
print(df_without_controls.shape[0])

# get positions
def get_chr(f):
    oligo_split = re.split('_|:|-', f)    
    chromosome = oligo_split[2]
    return chromosome
def get_start(f):
    oligo_split = re.split('_|:|-', f)    
    start = oligo_split[3]
    return start
def get_end(f):
    oligo_split = re.split('_|:|-', f)    
    end = oligo_split[4]
    return end
    
df_without_controls["chromosome"] = df_without_controls.oligo_column.map(get_chr)
df_without_controls["start+1"] = df_without_controls.oligo_column.map(get_start)
df_without_controls["end"] = df_without_controls.oligo_column.map(get_end)

df_without_controls["start"] = df_without_controls["start+1"].astype(int) -1


# remove duplicates due to L4
print(len(df_without_controls))
df_without_controls.drop(columns=['oligo_column'], inplace=True)
df_without_controls.drop_duplicates(keep = 'last', inplace=True)
print(len(df_without_controls))

print(df_without_controls.head().to_string())
df_without_controls.to_csv(f'./oligo_fasta/bed_files/L1L2L3L4a1_wo_controls.bed', header=True, index = False, columns=["chromosome", "start", "end"], sep='\t') #Added L4a1