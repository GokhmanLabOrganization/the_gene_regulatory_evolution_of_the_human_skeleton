# switch the replicates in barcode_counts_UMI.txt

import sys
import pandas as pd

# library = sys.argv[1]
# adaptor = sys.argv[2]
# cells = sys.argv[3]

# umi_file = f"/home/labs/davidgo/Collaboration/humanMPRA/{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI.txt"
umi_file = f"/home/labs/davidgo/Collaboration/RepCarlyClean/output/UMI/barcode_counts_UMI.txt"


umi_df = pd.read_csv(umi_file, sep='\t', index_col=0)

# umi_df[['RNA_chondrocytes_rep1','RNA_chondrocytes_rep2','DNA_chondrocytes_rep1','DNA_chondrocytes_rep2']]=umi_df[['RNA_chondrocytes_rep2','RNA_chondrocytes_rep1','DNA_chondrocytes_rep2','DNA_chondrocytes_rep1']]
umi_df[['RNA_Hob_rep1','RNA_Hob_rep2','DNA_Hob_rep1','DNA_Hob_rep2']]=umi_df[['RNA_Hob_rep2','RNA_Hob_rep1','DNA_Hob_rep2','DNA_Hob_rep1']]
umi_df_hob = umi_df.copy().filter(regex="Hob")

# umi_df.to_csv(f'/home/labs/davidgo/Collaboration/humanMPRA/additional/switching_replicates/{cells}/{library}{adaptor}/output/UMI/barcode_counts_UMI.txt', sep='\t', header=True)
umi_df_hob.to_csv(f'/home/labs/davidgo/Collaboration/RepCarlyClean/additional/switching_replicates_new_script/Hob/L0a0/output/UMI/barcode_counts_UMI_switched.txt', sep='\t', header=True)
