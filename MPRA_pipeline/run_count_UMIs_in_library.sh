#!/bin/bash
# bsub -q gsla-cpu -R "rusage[mem=40000]" -e ./additional/L1a1_UMI_counts_before_after/log/count_UMIs_inlibrary_%J.e.txt -o ./additional/L1a1_UMI_counts_before_after/log/count_UMIs_inlibrary_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_count_UMIs_in_library.sh $library $adaptor $cells

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/count_UMIs_in_library.py $1 $2 $3