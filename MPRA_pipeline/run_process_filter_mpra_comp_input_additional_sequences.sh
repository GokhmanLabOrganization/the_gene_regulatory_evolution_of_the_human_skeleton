#!/bin/bash
# cells=chondrocytes
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/comparative_analysis_combined/log/filter_mpra_comp_input_additional_sequences_%J.e.txt -o ./$cells/comparative_analysis_combined/log/filter_mpra_comp_input_additional_sequences_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_filter_mpra_comp_input_additional_sequences.sh $cells

module load python/3.7.9

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/filter_mpra_comp_input_additional_sequences.py $1