# this script filters the input of MPRAnalyze quantitative to only contain oligos that have more than 500 DNA counts

import pandas as pd
import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

RNA_org = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_quantitative/{cells}_{library}{adaptor}_RNA_MPRAnalyze_quantitative.txt', sep='\t', header =0)
DNA_org = pd.read_csv(f'./{cells}/{library}{adaptor}/output/mpranalyze_quantitative/{cells}_{library}{adaptor}_DNA_MPRAnalyze_quantitative.txt', sep='\t', header =0)

mask = DNA_org.sum(axis=1) >= 500

DNA_filtered = DNA_org[mask].copy()

RNA_filtered = RNA_org[mask].copy()

print(DNA_filtered.shape[0])

RNA_filtered.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_quantitative/{cells}_{library}{adaptor}_RNA_MPRAnalyze_quantitative_filter_DNA.txt', sep='\t', header=True, index=False)
DNA_filtered.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/mpranalyze_quantitative/{cells}_{library}{adaptor}_DNA_MPRAnalyze_quantitative_filter_DNA.txt', sep='\t', header=True, index=False)