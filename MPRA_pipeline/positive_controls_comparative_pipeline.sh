cells=chondrocytes #Fibroblastindels
library=L3
adaptor=a2
min_associations=2


cd /home/labs/davidgo/Collaboration/humanMPRA/

# MPRAnalyze comparative
# refer to MPRAnalyze comparative.docx


project=posctrl
mkdir -p ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/


bsub -q gsla-cpu -R "rusage[mem=6400]" -e ./$cells/$library$adaptor/log/format_for_MPRAnalyze_comparative_analysis_posctrl_%J.e.txt -o ./$cells/$library$adaptor/log/format_for_MPRAnalyze_comparative_analysis_posctrl_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_format_for_MPRAnalyze_comparative_analysis_posctrl.sh $library $adaptor $cells $project

#filter (runtime - a few minutes) - Look at the filtered counts in log file. 
bsub -q gsla-cpu -R "rusage[mem=6400]" -e ./$cells/$library$adaptor/log/filter_mpra_comp_input_for_active_posctrl_%J.e.txt -o ./$cells/$library$adaptor/log/filter_mpra_comp_input_for_active_posctrl_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_filter_mpra_comp_input_for_active_posctrl.sh $library $adaptor $cells $project 

#sort (runtime - 0.5 minute) NOTE: Make sure in script that additional_sequences is False!!
bsub -q gsla-cpu -R "rusage[mem=1600]" -e ./$cells/$library$adaptor/log/sort_comp_input_barcodes_posctrl_%J.e.txt -o ./$cells/$library$adaptor/log/sort_comp_input_barcodes_posctrl_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_sort_comp_input_barcodes_posctrl.sh $library $adaptor $cells $project

#csplit (runtime - ~0 minutes) NOTE: Make sure in script that additional_sequences is False!!
bsub -q gsla-cpu -R "rusage[mem=1600]" -e ./$cells/$library$adaptor/log/comp_csplit_posctrl_%J.e.txt -o ./$cells/$library$adaptor/log/comp_csplit_posctrl_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/csplit_posctrl.sh $library $adaptor $cells $project
input=filter_sorted
mkdir -p ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/res_filter_sorted/

for f in /home/labs/davidgo/Collaboration/humanMPRA/${cells}/${library}${adaptor}/output/mpranalyze_comparative/${project}/dna_chunks_${input}/*; do dnaname=$(basename $f ); rnaname="rna"${dnaname:3}; /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_for_loop_posctrl.sh $dnaname $rnaname ${dnaname:3} $library $adaptor $cells $input $project; done

cd  ./$cells/$library$adaptor/output/mpranalyze_comparative/$project/

#Join chunks directly from shell. Make sure that no additional seq oligos are also in that folder!!!
head -n 1 res_filter_sorted/001.txt > mpranalyze_comp_res_filter_sorted.txt
tail -n +2 -q res_filter_sorted/*.txt >> mpranalyze_comp_res_filter_sorted.txt
