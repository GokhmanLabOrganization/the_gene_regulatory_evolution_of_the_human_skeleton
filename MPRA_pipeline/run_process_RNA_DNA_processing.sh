#!/bin/bash
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/run_process_RNA_DNA_processing_%J.e.txt -o ./$cells/$library$adaptor/log/run_process_RNA_DNA_processing_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_RNA_DNA_processing.sh $library $adaptor $cells $min_associations

#module load python/3.7.9
pip install pandas --user
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/RNA_DNA_processing.py $1 $2 $3 $4
