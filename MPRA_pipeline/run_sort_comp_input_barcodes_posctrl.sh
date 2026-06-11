#!/bin/bash
# cells=chondrocytes
# library=L3
# adaptor=a1


# for library in L1 L2 L3;
# do for adaptor in a1 a2 a3;
# do cells=chondrocytes
# bsub -q new-short -R "rusage[mem=40000]" -e ./$cells/$library$adaptor/log/sort_comp_input_barcodes_%J.e.txt -o ./$cells/$library$adaptor/log/sort_comp_input_barcodes_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_sort_comp_input_barcodes.sh $library $adaptor $cells
# done
# done

#module load python/3.7.9

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/sort_comp_input_barcodes_posctrl.py $1 $2 $3 $4