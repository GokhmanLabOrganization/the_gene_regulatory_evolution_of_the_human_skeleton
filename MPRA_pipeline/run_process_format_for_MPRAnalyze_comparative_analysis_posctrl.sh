#!/bin/bash
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/format_for_MPRAnalyze_comparative_analysis_%J.e.txt -o ./$cells/$library$adaptor/log/format_for_MPRAnalyze_comparative_analysis_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_format_for_MPRAnalyze_comparative_analysis.sh $library $adaptor $cells $min_associations
module load python/3.7.9

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/format_for_MPRAnalyze_comparative_analysis_final_posctrl.py $1 $2 $3 $4
