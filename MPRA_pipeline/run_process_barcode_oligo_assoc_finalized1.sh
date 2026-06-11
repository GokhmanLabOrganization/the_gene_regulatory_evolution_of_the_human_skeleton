#!/bin/bash
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/run_process_barcode_oligo_assoc_finalized1_%J.e.txt -o ./$cells/$library$adaptor/log/run_process_barcode_oligo_assoc_finalized1_%J.o.txt ./run_process_barcode_oligo_assoc_finalized1.sh
module load python/3.7.9
pip install pandas --user
pip install pysam --user

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/process_barcode_oligo_assoc_finalized1.py $1 $2 $3
