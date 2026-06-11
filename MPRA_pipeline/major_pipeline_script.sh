# cells=neurons
# library=L3
# adaptor=a3
# min_associations=2


cd /home/labs/davidgo/Collaboration/humanMPRA/

################################## DNA barcode ##################################

bsub -q gsla-cpu -e ./$cells/$library$adaptor/log/map_barcodes_comb_finalized_%J.e.txt -o ./$cells/$library$adaptor/log/map_barcodes_comb_finalized1_%J.o.txt -n 16 -R "rusage[mem=8000]" -R "span[hosts=1]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/map_barcodes_comb_finalized1.sh $library $adaptor $cells

bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/run_process_barcode_oligo_assoc_finalized1_%J.e.txt -o ./$cells/$library$adaptor/log/run_process_barcode_oligo_assoc_finalized1_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_barcode_oligo_assoc_finalized1.sh $library $adaptor $cells

bsub -q gsla-cpu -R "rusage[mem=8000]" -e ./$cells/$library$adaptor/log/run_txt_to_fasta_finalized1_%J.e.txt -o ./$cells/$library$adaptor/log/run_txt_to_fasta_finalized1_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_txt_to_fasta_finalized1.sh $library $adaptor $cells $min_associations


########### DNA RNA (Run time: bowtie - ~4h, processing ~0.5h) ##################
bsub -q gsla-cpu -e ./$cells/$library$adaptor/log/alignBowtie_RNADNA_loop_%J.e.txt -o ./$cells/$library$adaptor/log/alignBowtie_RNADNA_loop_%J.o.txt -n 16 -R "rusage[mem=8000]" -R "span[hosts=1]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/alignBowtie_RNADNA_loop.sh $library $adaptor $cells $min_associations

bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/run_process_RNA_DNA_processing_%J.e.txt -o ./$cells/$library$adaptor/log/run_process_RNA_DNA_processing_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_RNA_DNA_processing.sh $library $adaptor $cells $min_associations

# filter for ouliers - optional (~5 minutes)
bsub -q short  -R "rusage[mem=50000]" -e ./$cells/$library$adaptor/log/filter_outliers_std2_%J.e.txt -o ./$cells/$library$adaptor/log/filter_outliers_std2_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_filter_outliers_std2.sh $library $adaptor $cells

# format for MPRAnalyze (run time quant - ~30 min, run time comp ~1hour) 

bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/format_for_MPRAnalyze_quantitative_analysis_%J.e.txt -o ./$cells/$library$adaptor/log/format_for_MPRAnalyze_quantitative_analysis_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_format_for_MPRAnalyze_quantitative_analysis.sh $library $adaptor $cells

bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/format_for_MPRAnalyze_comparative_analysis_%J.e.txt -o ./$cells/$library$adaptor/log/format_for_MPRAnalyze_comparative_analysis_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_format_for_MPRAnalyze_comparative_analysis.sh $library $adaptor $cells

################################## MPRAnalyze Quant ##################################

# MPRAnalyze quantitative (~ 1hour)

bsub -q gsla-cpu -R "rusage[mem=200000]" -e ./$cells/$library$adaptor/log/run_process_MPRAnalyze_quantitative_%J.e.txt -o ./$cells/$library$adaptor/log/run_process_MPRAnalyze_quantitative_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_quantitative.sh $library $adaptor $cells

# calculate ratio

bsub -q gsla-cpu -R "rusage[mem=50000]" -e ./$cells/$library$adaptor/log/run_calculate_ratio_df_%J.e.txt -o ./$cells/$library$adaptor/log/run_calculate_ratio_df_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_calculate_ratio_df.sh $library $adaptor $cells $min_associations

# quant results (runtime - about 20 min)

bsub -q gsla-cpu -R "rusage[mem=40000]" -e ./$cells/$library$adaptor/log/run_quant_results_%J.e.txt -o ./$cells/$library$adaptor/log/run_quant_results_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_quant_results.sh $library $adaptor $cells $min_associations

################################## MPRAnalyze comp ##################################
# refer to MPRAnalyze comparative.docx

#filter (runtime - a few minutes) - Look at the filtered counts in log file. 
bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/filter_mpra_comp_input_for_active_%J.e.txt -o ./$cells/$library$adaptor/log/filter_mpra_comp_input_for_active_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_filter_mpra_comp_input_for_active.sh $library $adaptor $cells

#sort (runtime - 0.5 minute) NOTE: Make sure in script that additional_sequences is False!!
bsub -q gsla-cpu -R "rusage[mem=40000]" -e ./$cells/$library$adaptor/log/sort_comp_input_barcodes_%J.e.txt -o ./$cells/$library$adaptor/log/sort_comp_input_barcodes_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_sort_comp_input_barcodes.sh $library $adaptor $cells

#csplit (runtime - ~0 minutes) NOTE: Make sure in script that additional_sequences is False!!
bsub -q gsla-cpu -R "rusage[mem=40000]" -e ./$cells/$library$adaptor/log/comp_csplit_%J.e.txt -o ./$cells/$library$adaptor/log/comp_csplit_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/csplit.sh $library $adaptor $cells

#MPRAnalyze loop
# cells=chondrocytes
# library=L3
# adaptor=a3
# input=filter_sorted

#NOTE: Make sure in script that additional_sequences is False!!
for f in /home/labs/davidgo/Collaboration/humanMPRA/${cells}/${library}${adaptor}/output/mpranalyze_comparative/dna_chunks_${input}/*; do dnaname=$(basename $f ); rnaname="rna"${dnaname:3}; /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_for_loop.sh $dnaname $rnaname ${dnaname:3} $library $adaptor $cells $input; done

#Join chunks directly from shell. Make sure that no additional seq oligos are also in that folder!!!
head -n 1 res_filter_sorted/001.txt > mpranalyze_comp_res_filter_sorted.txt
tail -n +2 -q res_filter_sorted/*.txt >> mpranalyze_comp_res_filter_sorted.txt #Make sure that no additional seq oligos are also in that folder!!!


##################################Combined analysis##################################

# 1. Combine quant results of all relevant libries for the specific cell type (about 1.5 minutes of runtime)

#cells=chondrocytes
bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./$cells/quantitative_analysis_combined/log/combining_quant_results_%J.e.txt -o ./$cells/quantitative_analysis_combined/log/combining_quant_results_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_combining_quant_results.sh $cells

# 2. combine comparative results. NOTE: make sure include_additional_sequences is False!
bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./$cells/comparative_analysis_combined/log/add_fdr_to_mpranalyze_comp_%J.e.txt -o ./$cells/comparative_analysis_combined/log/add_fdr_to_mpranalyze_comp_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_add_fdr_to_mpranalyze_comp.sh $cells

# 3. prepare comaprative for additional sequences after combinging FDR
bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/comparative_analysis_combined/log/filter_mpra_comp_input_additional_sequences_%J.e.txt -o ./$cells/comparative_analysis_combined/log/filter_mpra_comp_input_additional_sequences_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_filter_mpra_comp_input_additional_sequences.sh $cells

# 4. again run the sort, csplit, loop comparative with additional_sequences set to true
#csplit...

for f in /home/labs/davidgo/Collaboration/humanMPRA/${cells}/${library}${adaptor}/output/mpranalyze_comparative/dna_chunks_${input}/dna_additional_sequences*; do dnaname=$(basename $f ); rnaname="rna"${dnaname:3}; /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_for_loop.sh $dnaname $rnaname ${dnaname:3} $library $adaptor $cells $input; done

# 5. For each library, combined the new additional sequences into a 
head -n 1 res_filter_sorted/_additional_sequences001.txt > additional_sequences_mpranalyze_comp_res_filter_sorted.txt
tail -n +2 -q res_filter_sorted/_additional_sequences*.txt >> additional_sequences_mpranalyze_comp_res_filter_sorted.txt

# 6. combine comparative results. NOTE: make sure include_additional_sequences is True! (runtime = ~2m)
bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./$cells/comparative_analysis_combined/log/add_fdr_to_mpranalyze_comp_%J.e.txt -o ./$cells/comparative_analysis_combined/log/add_fdr_to_mpranalyze_comp_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_add_fdr_to_mpranalyze_comp.sh $cells


# 7.combine final quantitative and final comparative tables (runtime ~1 min)

bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./$cells/comparative_analysis_combined/log/combine_quant_comp_%J.e.txt -o ./$cells/comparative_analysis_combined/log/combine_quant_comp_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_combine_quant_comp.sh $cells