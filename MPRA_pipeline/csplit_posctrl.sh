#!/bin/bash

library=$1
adaptor=$2
cells=$3  
project=$4

mkdir -p ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/rna_chunks_filter_sorted/
mkdir -p ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/dna_chunks_filter_sorted/

csplit -f ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/rna_chunks_filter_sorted/rna -s -n 3 ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/${cells}_${library}${adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt $(cat ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/split_positions.txt)
csplit -f ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/dna_chunks_filter_sorted/dna -s -n 3 ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/${cells}_${library}${adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted.txt $(cat ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/split_positions.txt)


## for additional sequences: 

#csplit -f ./$cells/$library$adaptor/output/mpranalyze_comparative/rna_chunks_filter_sorted/rna_additional_sequences -s -n 3 ./$cells/$library$adaptor/output/mpranalyze_comparative/${cells}_${library}${adaptor}_RNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted_additional_sequences.txt $(cat ./$cells/$library$adaptor/output/mpranalyze_comparative/split_positions_additional_sequences.txt)
#csplit -f ./$cells/$library$adaptor/output/mpranalyze_comparative/dna_chunks_filter_sorted/dna_additional_sequences -s -n 3 ./$cells/$library$adaptor/output/mpranalyze_comparative/${cells}_${library}${adaptor}_DNA_MPRAnalyze_comparative_filter_adjusted_fdr_sorted_additional_sequences.txt $(cat ./$cells/$library$adaptor/output/mpranalyze_comparative/split_positions_additional_sequences.txt)
