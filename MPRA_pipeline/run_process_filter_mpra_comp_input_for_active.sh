#!/bin/bash
# cells=chondrocytes
# library=L3
# adaptor=a1


# for library in L1 L2 L3;
# do adaptor=a3;
# cells=chondrocytes
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/filter_mpra_comp_input_for_active_%J.e.txt -o ./$cells/$library$adaptor/log/filter_mpra_comp_input_for_active_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_filter_mpra_comp_input_for_active.sh $library $adaptor $cells
# done

#module load python/3.7.9

python3 /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/filter_mpra_comp_input_for_active.py $1 $2 $3