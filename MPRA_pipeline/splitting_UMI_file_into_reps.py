# split merged UMI file into replicates
import pandas as pd
import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

UMI_df = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI.txt', sep='\t', header=0)
for rep in ["rep1", "rep2", "rep3"]:
	UMI_df.to_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI_{rep}.txt', sep='\t', header=True, index=False, columns=['oligo_bc', f'RNA_{cells}_{rep}', f'DNA_{cells}_{rep}'])