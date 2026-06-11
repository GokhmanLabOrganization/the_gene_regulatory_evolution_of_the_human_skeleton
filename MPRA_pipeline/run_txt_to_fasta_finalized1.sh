#!/bin/bash
# bsub -q molgen-q -R "rusage[mem=8000]" -e ./$cells/$library$adaptor/log/run_txt_to_fasta_finalized1_%J.e.txt -o ./$cells/$library$adaptor/log/run_txt_to_fasta_finalized1_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_txt_to_fasta_finalized1.sh $library $adaptor $cells $min_associations
module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/txt_to_fasta_finalized1.py $1 $2 $3 $4
