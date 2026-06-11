#!/usr/bin/env python
#
# need to try to have an indexed bam with the input (prepared upstream) # Simon 4.7.21
import pandas as pd
import numpy as np
import sys
import pysam
import os
import re
import time
import cProfile
import matplotlib # added to prevent display error #Katharina 28.7.22
matplotlib.use('Agg') # added to prevent display error #Katharina 28.7.22
import matplotlib.pyplot as plt
from collections import Counter

library = sys.argv[1]
adaptor = sys.argv[2]
cells = sys.argv[3]

minimum_barcode_sequence_quality = 30
minimum_number_reads_associating_barcode_to_oligo = 2
minimum_oligo_mapping_quality = 6

print(minimum_barcode_sequence_quality)
print(minimum_number_reads_associating_barcode_to_oligo)

#project_dir = f'./{cells}/{library}{adaptor}/output/DNA_barcode_associations_2/'
project_dir = f'./{cells}/{library}{adaptor}/output/DNA_barcode_associations_2/Omer_re-send/' #temp for L4a1
fastq_dir = f'./{cells}/{library}{adaptor}/input/DNA_barcode'
fastq_file = f'{library}-{adaptor}_R2.fastq.gz'
bam_filename = f'./{cells}/{library}{adaptor}/output/Aligned_DNA_barcodes_comb.bam'
filtered_bam_filename = f'Aligned_DNA_barcodes_comb_filtered_{library}{adaptor}.bam'
barcode_txt_doc = f'oligos_to_barcodes_comb_{library}{adaptor}_{minimum_barcode_sequence_quality}_{minimum_number_reads_associating_barcode_to_oligo}.txt'

def turn_bam_into_molecules(bam_file):
    
    ## TO DO ##
    # filter for proper_pair (if false for read1, toss both reads)
    # also filter for mapq >=6, check for read 1 and read 2 separately, keep either if 6 or more
    # if we want to be more strict, we can check quality of variant base at either or both (whichever passed so far), only keep reads where quality of that base is equal to or more than 30 (we don't want to filter on "them being associated with the correct oligo" 3 times)
    
    bam_molecule_dict_filtered = {}
    not_proper_pair = 0
    molecule_names_didnt_match = 0
    not_mapping_to_same_oligo = 0
    number_molecules_in_filtered_bam = 0
    number_read1_bad_map = 0
    number_read2_bad_map = 0
    
    # copied/inspired from JJ script map_barcodes.py
    
    with pysam.AlignmentFile(f'{project_dir}/{filtered_bam_filename}', 'wb', template = bam_file) as outfile:
    
        for i, read1 in enumerate(bam_file):

            if not read1.is_read1:  # go through reads as pairs only, skip every other read (don't look at read 2 if read 1 not present)
                continue
            
            if not read1.is_proper_pair:  # filter out cases that aren't 'proper pairs' (make sure paired-endedness worked ok?)
                not_proper_pair += 1
                continue
            
            read2 = next(bam_file)
            if not read2.is_read2:
                continue # make sure that we're operating on pairs of reads
            
            # make sure the molecule names match each other
            read1_molecule = read1.query_name
            read2_molecule = read2.query_name
            
            if not read1_molecule == read2_molecule:  
                molecule_names_didnt_match += 1
                continue
            
            # make sure they map to the same oligo
            read1_oligo = read1.reference_name
            read2_oligo = read2.reference_name
            
            if not read1_oligo == read2_oligo:  
                not_mapping_to_same_oligo += 1
                continue
            
            # filter on mapq and write to new filtered bam
            read1_mapping_quality = read1.mapping_quality
            read2_mapping_quality = read2.mapping_quality
            
            if read1_mapping_quality >= minimum_oligo_mapping_quality:
                number_molecules_in_filtered_bam += 1
                outfile.write(read1)
            if read1_mapping_quality < minimum_oligo_mapping_quality:
                number_read1_bad_map += 1
            
            if read2_mapping_quality >= minimum_oligo_mapping_quality:
                number_molecules_in_filtered_bam += 1
                outfile.write(read2)
            if read2_mapping_quality < minimum_oligo_mapping_quality:
                number_read2_bad_map += 1
            
            # store the molecule in a dict as long as one of them passes mapq filtering
            if read1_mapping_quality >= minimum_oligo_mapping_quality or read2_mapping_quality >= minimum_oligo_mapping_quality:
                bam_molecule_dict_filtered[read1_molecule] = [read1_oligo]
   
    number_original_alignments = pysam.view("-c", bam_filename)
    number_filtered_alignments = pysam.view("-c", f'{project_dir}/{filtered_bam_filename}')
    
    print(number_original_alignments, 'number of original alignments')
    print(number_filtered_alignments, 'number of filtered alignments')
    print(number_molecules_in_filtered_bam, 'number of alignments in bam after all this filtering, should be same as above')
    print(not_proper_pair, 'number of pairs of alignments filtered due to not proper pair')
    print(molecule_names_didnt_match, 'number of alignment pairs filtered because molecules didnt match between reads')
    print(not_mapping_to_same_oligo, 'number of alignment pairs filtered because the two reads didnt map to same oligo')
    print(len(bam_molecule_dict_filtered),'number of unique molecules after filtering (includes R1 and R3)')
    print(number_read1_bad_map, 'number read1 filtered for bad mapping')
    print(number_read2_bad_map, 'number read2 filtered for bad mapping')

    return bam_molecule_dict_filtered   # each molecule (key) has a single value, the oligo it mapped to. molecules will not be repeated, but oligos will be 

def turn_fastq_into_barcode_dict(fastq, filtered_bam_molecule_dict):

    fastq_barcode_dict={}
    total_entries = 0
    number_entries_not_found_in_filtered_bam = 0
    number_barcodes_not_good_enough_quality = 0
    entries_going_into_dict = 0

    for each_entry in fastq:
        total_entries += 1
        
        # only use molecules that are in the filtered bam (we won't have oligo info for the molecules that aren't anyways)
        molecule = each_entry.name
        if not molecule in filtered_bam_molecule_dict:
            number_entries_not_found_in_filtered_bam += 1
            continue
        
        entries_going_into_dict += 1  # count up how many times we are updating the bam molecule dict
        
        barcode_sequence = each_entry.sequence
        min_base_quality = min(each_entry.get_quality_array())  # ints, already have 33 subtracted
        
        if min_base_quality < minimum_barcode_sequence_quality:   # if we want to filter out the barcodes with a certain minimum base quality we can add continue here
            number_barcodes_not_good_enough_quality += 1
        
        filtered_bam_molecule_dict[molecule].extend([barcode_sequence, min_base_quality])  # add info to the original bam dict. now molecule: [oligo, barcode_sequence, min_base_quality]
        
        if barcode_sequence in fastq_barcode_dict:
            fastq_barcode_dict[barcode_sequence].append([molecule, min_base_quality])   # ends up as barcode: [[mol1, min_base_qual1], [mol2, min_base_qual2]]
        else:
            fastq_barcode_dict[barcode_sequence] = [[molecule, min_base_quality]]
    
    mol_count = sum(len(v) for v in fastq_barcode_dict.values())
    unique_bc = len(fastq_barcode_dict)
    
    print(total_entries, 'total entries in fastq')
    print(number_entries_not_found_in_filtered_bam, 'number of fastq entries not found in filtered bam')
    print(number_barcodes_not_good_enough_quality, 'number of barcodes with minimum base quality not good enough (but not filtered out)')
    print(mol_count, 'total molecules in final barcode_dict')
    print(unique_bc, 'total unique barcodes in barcode dict, not qual filtered for min base qual')
    print(entries_going_into_dict, 'fastq entries that updated the bam dict')
    
    return fastq_barcode_dict, filtered_bam_molecule_dict    # for each barcode (key), have list of all the molecules in the fastq that mapped that same barcode at high quality. will be many. potentially have molecules mapping to different oligos, will filter later

def turn_dict_into_dataframe(everything_dict):
    
    dataframe = pd.DataFrame.from_dict(everything_dict, orient='index', columns=['oligo', 'bc_seq', 'bc_min_qual'])
    print(len(dataframe.index), 'number of mols in original dataframe')
    unique_bc_start = len(list(dataframe.bc_seq.unique()))
    print(unique_bc_start, 'number of unique barcodes before bcseq qual filter')
    
    #dataframe.to_csv(f'{project_dir}/df_mapq_filtered.txt', header=True)
    
    qual_hist = dataframe.hist(column='bc_min_qual', bins=20)
    plt.savefig(f'{project_dir}/bc_quality_hist_comb_{library}{adaptor}.pdf')
    
    
    # also filter for bc min qual here and then remove that column from df
    
    filtered_df = dataframe.loc[dataframe['bc_min_qual'] >= minimum_barcode_sequence_quality]
    smaller_df = filtered_df.drop(['bc_min_qual'], axis=1)
    unique_bc_fil1 = len(list(smaller_df.bc_seq.unique()))
    
    print(len(smaller_df.index), 'number of mols in bc qual filtered dataframe')
    print(unique_bc_fil1, 'number of unique barcodes after bcseq qual filter')
    return smaller_df

def filter_promisc_barcodes(df):
    
    barcodes_not_mapping_to_multiple_oligos = []
    promisc_barcodes = 0
    number_bc_with_2_or_more = 0
    
    # find barcodes that don't map to multiple oligos and only keep them
    df_bc_grouped = df.groupby('bc_seq')
    
    for each_barcode, barcode_group in df_bc_grouped:
        number_associations_seen = len(barcode_group.index)
        if number_associations_seen >= 2:
            number_bc_with_2_or_more += 1
        oligo_list = barcode_group['oligo'].to_list()
        if all(x==oligo_list[0] for x in oligo_list):
            barcodes_not_mapping_to_multiple_oligos.append(each_barcode)
        else:
            promisc_barcodes += 1
    
    
    print(len(barcodes_not_mapping_to_multiple_oligos), 'barcodes not mapping to multiple oligos, before minimum associations filter')
    print(promisc_barcodes, 'number of unique promiscuous barcodes')
    print(len(df), 'length of df before filtering for multiple oligos')
    print(number_bc_with_2_or_more, 'number of barcodes with at least 2 mapping to them before promisc. filter')
    filtered_df = df[df['bc_seq'].isin(barcodes_not_mapping_to_multiple_oligos)]  # only keeps rows where the value of bc_seq is in the list
    print(len(filtered_df), 'length of df after filtering for multiple oligos')
    
    return filtered_df
    
def filter_for_minimum_barcode_associations(df):
    
    barcodes_with_enough_associations = []
    number_reads_per_barcode = []
    
    df_bc_grouped = df.groupby('bc_seq')
    
    for each_barcode, barcode_group in df_bc_grouped:
        number_associations_seen = len(barcode_group.index)
        number_reads_per_barcode.append(number_associations_seen)
        if number_associations_seen >= minimum_number_reads_associating_barcode_to_oligo:
            barcodes_with_enough_associations.append(each_barcode)  # save good barcodes with enough associations
    
    plt.clf()
    plt.hist(number_reads_per_barcode)
    plt.savefig(f'{project_dir}/number_reads_per_barcode_hist_test_{library}{adaptor}.pdf')
    plt.clf()
    plt.hist(number_reads_per_barcode, bins=[0,1,2,3,4,5,6,7,8,9,10,15,20,50])
    plt.savefig(f'{project_dir}/number_reads_per_barcode_hist_{library}{adaptor}.pdf')
    plt.clf()
    plt.hist(number_reads_per_barcode, bins=[0,1,2,3,4,5])
    plt.savefig(f'{project_dir}/number_reads_per_barcode_hist_small_{library}{adaptor}.pdf')
    
    print(len(df), 'length of df before filtering for barcode associations')
    print(len(barcodes_with_enough_associations), 'number of final unique barcodes')
    filtered_df = df[df['bc_seq'].isin(barcodes_with_enough_associations)] # only keep rows where the value of bc_seq is in the list of barcodes with enough associations
    print(len(filtered_df), 'length of df after filtering for minimum barcode associations')
    #unprint when done:
    #filtered_df.to_csv(f'{project_dir}/final_filtered_barcode_reads_{minimum_number_reads_associating_barcode_to_oligo}_{library}{adaptor}.txt', index=None, header=True, sep='\t')
    return filtered_df
 
def save_barcode_oligo_mappings(df):
    
    unique_bc_oligo_mappings = df.drop_duplicates()
    unique_oligos = len(list(df.oligo.unique()))
    print(unique_oligos, 'number of unique oligos in final bc list')
    #unprint when done:
    #unique_bc_oligo_mappings.to_csv(f'{project_dir}/{barcode_txt_doc}', index=None, header=True)
    
def turn_bam_into_molecules_no_filter(bam_file):
    #This function was copied from Omer's code on 03-08-2025
    bam_molecule_dict = {}
    not_proper_pair = 0
    molecule_names_didnt_match = 0
    not_mapping_to_same_oligo = 0
    number_molecules_in_filtered_bam = 0
    number_read1_bad_align = 0
    number_read2_bad_align = 0

    # copied/inspired from JJ script map_barcodes.py

    with pysam.AlignmentFile(f'{project_dir}/{filtered_bam_filename}', 'wb', template=bam_file) as outfile:

        for i, read1 in enumerate(bam_file):

            if not read1.is_read1:  # go through reads as pairs only, skip every other read (don't look at read 2 if read 1 not present)
                continue

            if not read1.is_proper_pair:  # filter out cases that aren't 'proper pairs' (make sure paired-endedness worked ok?)
                not_proper_pair += 1
                continue

            read2 = next(bam_file)
            if not read2.is_read2:
                continue  # make sure that we're operating on pairs of reads

            # make sure the molecule names match each other
            read1_molecule = read1.query_name
            read2_molecule = read2.query_name

            if not read1_molecule == read2_molecule:
                molecule_names_didnt_match += 1
                continue

            # make sure they map to the same oligo
            read1_oligo = read1.reference_name
            read2_oligo = read2.reference_name

            if not read1_oligo == read2_oligo:
                not_mapping_to_same_oligo += 1
                continue

            # filter on mapq and write to new filtered bam
            read1_cigar = read1.cigarstring
            read2_cigar = read2.cigarstring
            read1_md = read1.get_tag('MD')
            read2_md = read2.get_tag('MD')

            if read1_cigar == "146M" and read1_md == '146':
                number_molecules_in_filtered_bam += 1
                outfile.write(read1)
            else:
                number_read1_bad_align += 1

            if read2_cigar == "146M" and read2_md =='146':
                number_molecules_in_filtered_bam += 1
                outfile.write(read2)
            else:
                number_read2_bad_align += 1

            # store the molecule in a dict as long as one of them passes mapq filtering
            #if (read1_cigar == "146M" and read1_md == '146') or (read2_cigar == "146M" and read2_md =='146'):

            bam_molecule_dict[read1_molecule] = [read1_oligo]

    number_original_alignments = pysam.view("-c", bam_filename)
    number_filtered_alignments = pysam.view("-c", f'{project_dir}/{filtered_bam_filename}')

    print(number_original_alignments, 'number of original alignments')
    print(number_filtered_alignments, 'number of filtered alignments')
    print(number_molecules_in_filtered_bam,
          'number of alignments in bam after all this filtering, should be same as above')
    print(not_proper_pair, 'number of pairs of alignments filtered due to not proper pair')
    print(molecule_names_didnt_match, 'number of alignment pairs filtered because molecules didnt match between reads')
    print(not_mapping_to_same_oligo, 'number of alignment pairs filtered because the two reads didnt map to same oligo')
    print(len(bam_molecule_dict), 'number of unique molecules after filtering (includes R1 and R3)')
    print(number_read1_bad_align, 'number read1 filtered for bad alignment')
    print(number_read2_bad_align, 'number read2 filtered for bad alignment')

    return bam_molecule_dict  # each molecule (key) has a single value, the oligo it mapped to. molecules will not be repeated, but oligos will be

### Run the script ###
def main():
    print("step 0")
    bam = pysam.AlignmentFile(bam_filename, 'rb')
    full_mol_dict=turn_bam_into_molecules_no_filter(bam)
    fastq_object = pysam.FastxFile(f'{fastq_dir}/{fastq_file}')
    barcode_dict, better_mol_dict = turn_fastq_into_barcode_dict(fastq_object, full_mol_dict)
    before_filtering_df = pd.DataFrame.from_dict(better_mol_dict, orient='index',
                                                 columns=['oligo', 'bc_seq', 'bc_min_qual'])
    before_filtering_df.to_csv(f'{project_dir}/before_filtering.csv')
    bam = pysam.AlignmentFile(bam_filename, 'rb')

    print("step 1")
    molecule_dict = turn_bam_into_molecules(bam)
    print("step 2")
    fastq_object = pysam.FastxFile(f'{fastq_dir}/{fastq_file}')
    barcode_dict, better_mol_dict = turn_fastq_into_barcode_dict(fastq_object, molecule_dict)
    after_quality_df = pd.DataFrame.from_dict(better_mol_dict, orient='index', columns=['oligo', 'bc_seq', 'bc_min_qual'])
    after_quality_df.to_csv(f'{project_dir}/after_quality.csv')
    print("step 3")    
    
    filtered_df = turn_dict_into_dataframe(better_mol_dict)
    filtered_df.to_csv(f'{project_dir}/after_barcode_filter.csv')
    print("step 4")
    no_multiple_oligos = filter_promisc_barcodes(filtered_df)
    no_multiple_oligos.to_csv(f'{project_dir}/after_promiscuity_filter.csv')
    print("step 5")
    enough_associations = filter_for_minimum_barcode_associations(no_multiple_oligos)
    enough_associations.to_csv(f'{project_dir}/after_associations_filter.csv')
    print("step 6")
    save_barcode_oligo_mappings(enough_associations)
    print("step 7")
    
main()
    

