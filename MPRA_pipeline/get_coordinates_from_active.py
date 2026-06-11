# create bed file from active oligos

import pandas as pd
import re
import sys     
      
cells = sys.argv[1]

quant = pd.read_csv(f'./{cells}/quantitative_analysis_combined/comb_df_combined_fdr.csv', header=0)
print(quant.shape[0])
#filter for active
quant_active = quant[quant["activity_adjusted_combined"]=="active"] 
print(quant_active.shape[0])

# remove allele from oligo name
quant_active['oligo'] = quant_active['oligo'].str.replace("_ancestral|_derived", '')

# remove duplicates based on oligo column
quant_active.drop_duplicates(subset=['oligo'], inplace=True)
print(quant_active.shape[0])

# turn oligo column into bed file
def extract_and_format(column_value):
    # Split the string by underscores
    parts = column_value.split('_')
    # Extract the part between the second and third underscore
    extracted_part = parts[2]
    # Split the extracted part by ':' and '-'
    split_parts = re.split(':|-', extracted_part)
    # Create three separate columns and substract 1 from start to get 0 based
    col1, col2, col3 = split_parts[0], int(split_parts[1]) - 1, int(split_parts[2])

    return pd.Series([col1, col2, col3])

# Apply the function to the oligo column
coordinates = quant_active['oligo'].apply(extract_and_format)

# Add column names
coordinates.columns = ['chromosome', 'start', 'end']

coordinates.to_csv(f'./{cells}/quantitative_analysis_combined/active_sequences_coordinates.txt', sep='\t', index=False)