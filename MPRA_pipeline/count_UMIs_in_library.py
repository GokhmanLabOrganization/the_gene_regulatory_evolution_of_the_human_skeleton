import pandas as pd
import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

barcode_counts = pd.read_csv(f'./{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI.txt', sep='\t', header =0)
print(barcode_counts.head())

#Count the number of UMIs in rep1
RNA_counts = barcode_counts['RNA_chondrocytes_rep1'].sum()
print(RNA_counts)
DNA_counts = barcode_counts['DNA_chondrocytes_rep1'].sum()
print(DNA_counts)
total_counts = RNA_counts+DNA_counts
print(total_counts)

