import argparse
import pandas as pd
from pathlib import Path

# Load the parameters from the config file
parser = argparse.ArgumentParser(description="config file path")
parser.add_argument('-c', '--config_file_path', help='path to the config file')
args = parser.parse_args()
PARAMS = pd.read_csv(args.config_file_path, delimiter='\t', index_col=0, header=None).to_dict()[1]

# Function to merge duplicate IDs
def merge_gene_rows(group):
    merged_row = group.ffill().bfill().iloc[0]  # Fill missing values
    if group["chromosome"].isna().any():
        merged_row["Gene"] = group.loc[group["chromosome"].isna(), "Gene"].values[0]
    return merged_row

if __name__ == '__main__':
    ref_lfc = pd.read_csv(Path(PARAMS['ref_deseq2']), delimiter='\t')
    alt_lfc = pd.read_csv(Path(PARAMS['alt_deseq2']), delimiter='\t')

    ref_suffix = f"{PARAMS['ref_species']}_ref"
    alt_suffix = f"{PARAMS['alt_species']}_ref"

    combined_lfc = pd.merge(ref_lfc, 
                            alt_lfc, 
                            on="Gene", 
                            how="outer",
                            suffixes=(f"_{ref_suffix}", f"_{alt_suffix}"))
    
    # # Use NCBI gene annotations to fill in missing chromosome info (for aneuploidy) and GeneID info to merge synonyms symbols
    # NCBI_hg19 = pd.read_csv('/home/labs/davidgo/Collaboration/GenomeAnnotation/Human/NCBI/hg19/gff_parsed/NCBI_genes_hg19.txt', 
    #                     delimiter='\t')
    # combined_lfc = combined_lfc.merge(NCBI_hg19[['chromosome', 'GeneName', 'GeneID']], 
    #                                   left_on='Gene', 
    #                                   right_on='GeneName', 
    #                                   how='left', 
    #                                   suffixes=('', '_ref'))
    
    # merged_df = combined_lfc.groupby("GeneID", group_keys=False).apply(merge_gene_rows).reset_index(drop=True)
    # # Append rows where GeneID is NaN (no match found) back to the merged dataframe
    # combined_lfc = pd.concat([merged_df, combined_lfc[combined_lfc['GeneID'].isna()]])


    # Create masks for each filtering criterion
    na_mask = (pd.isna(combined_lfc[f'log2FoldChange_{ref_suffix}']) |
               pd.isna(combined_lfc[f'padj_{ref_suffix}']) |
               pd.isna(combined_lfc[f'log2FoldChange_{alt_suffix}']) |
               pd.isna(combined_lfc[f'padj_{alt_suffix}']))
    
    non_significant_mask = (combined_lfc[f'padj_{ref_suffix}'] > 0.05) | \
                           (combined_lfc[f'padj_{alt_suffix}'] > 0.05)
                     
    non_agreement_mask = (combined_lfc[f'log2FoldChange_{ref_suffix}'] *
                            combined_lfc[f'log2FoldChange_{alt_suffix}'] < 0)
                      
    non_close_value_mask = (combined_lfc[f'log2FoldChange_{ref_suffix}'] - 
                            combined_lfc[f'log2FoldChange_{alt_suffix}']).abs() >= 1

    # # Remove genes in aneuploid chromosomes if specified (do not need if they are removed prior to DESeq2 with TPM_calculation.py)
    # if PARAMS['aneuploidy'] != "":
    #     aneuploidy_chromosomes = [chrom.strip() for chrom in PARAMS['aneuploidy'].split(',')]
    #     aneuploidy_mask = combined_lfc['chromosome'].isin(aneuploidy_chromosomes)
    #     na_mask |= aneuploidy_mask

    # Calculate the final masks for each category
    rest_mask = na_mask | non_agreement_mask | non_close_value_mask
    non_ase_mask = (~rest_mask) & (non_significant_mask)
    ase_mask = ~(rest_mask | non_ase_mask)

    # Assign categories based on the masks on column 'ASE_definition'
    combined_lfc.loc[rest_mask, 'ASE_definition'] = 'else'
    combined_lfc.loc[non_ase_mask, 'ASE_definition'] = 'nonASE'
    combined_lfc.loc[ase_mask, 'ASE_definition'] = 'ASE'

    print(combined_lfc['ASE_definition'].value_counts())

    combined_lfc.to_csv(Path(PARAMS['output_path']) / 'combined_lfc_with_ASE_definition.tsv', sep='\t', index=False)