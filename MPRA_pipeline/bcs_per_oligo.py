# 08-10-2024 Nadav Mishol
# This script prints histograms of Barcodes per oligo for a specific library
# Original file by Katharina/Carly, Nadav copied it from Collaboration\backup\humanMPRA\scripts to the backup folder


import pandas as pd
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
import sys

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]
associations = sys.argv[4]

def hist_number_of_bcs_per_oligo(df):
    number_bcs_per_oligo = []
    df_oligo_grouped = df.groupby('oligo')
    for each_oligo, oligo_group in df_oligo_grouped:
        number_bcs_seen = len(oligo_group.index)
        number_bcs_per_oligo.append(number_bcs_seen)
        
        
    plt.clf()
    plt.hist(number_bcs_per_oligo)
    plt.savefig(f'./additional/bcs_per_oligo/number_bcs_per_oligo_hist_test_{cells}_{library}{adaptor}_{associations}.pdf')
    
    plt.clf()
    plt.hist(number_bcs_per_oligo, range = (0,500))
    plt.savefig(f'./additional/bcs_per_oligo/number_bcs_per_oligo_hist_0_500_{cells}_{library}{adaptor}_{associations}.pdf')
    
    plt.clf()
    plt.hist(number_bcs_per_oligo, range = (0,50))
    plt.savefig(f'./additional/bcs_per_oligo/number_bcs_per_oligo_hist_0_50_{cells}_{library}{adaptor}_{associations}.pdf')
    
    # plt.clf()
    # plt.hist(number_bcs_per_oligo, bins=[0,1,2,3,4,5,6,7,8,9,10,15,20,50,100])
    # plt.savefig(f'./additional/bcs_per_oligo/number_bcs_per_oligo_hist_{cells}_{library}{adaptor}_{associations}.pdf')
    
    count10 = len([elem for elem in number_bcs_per_oligo if elem >= 10])
    count5 = len([elem for elem in number_bcs_per_oligo if elem >= 5])
    count20 = len([elem for elem in number_bcs_per_oligo if elem >= 20])
    total = df_oligo_grouped.ngroups
    ratio10 = count10/total
    
    print('Count of oligos with at least 10 bcs: ', count10)
    print('Count of oligos with at least 5 bcs: ', count5)
    print('Count of oligos with at least 20 bcs: ', count20)
    print('Total oligos: ', total)
    print('fraction of oligos with at least 10 bcs: ', ratio10)
    print('max number of bcs per oligo: ', max(number_bcs_per_oligo))
    

oligo_bcs_df = pd.read_csv(f'./{cells}/{library}{adaptor}/output/DNA_barcode_associations_{associations}/oligos_to_barcodes_comb_{library}{adaptor}_30_{associations}.txt')
hist_number_of_bcs_per_oligo(oligo_bcs_df)