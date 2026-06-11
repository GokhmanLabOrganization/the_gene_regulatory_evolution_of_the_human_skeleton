#!/bin/bash

# cells=chondrocytes
# library=L3
# adaptor=a2
# input=seq_325030_v2

# for f in /home/labs/davidgo/Collaboration/humanMPRA/${cells}/${library}${adaptor}/output/mpranalyze_comparative/dna_chunks_${input}/*; do dnaname=$(basename $f ); rnaname="rna"${dnaname:3}; /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_for_loop.sh $dnaname $rnaname ${dnaname:3} $library $adaptor $cells $input; done

# # for additional sequences after correcting fdr:
# for f in /home/labs/davidgo/Collaboration/humanMPRA/${cells}/${library}${adaptor}/output/mpranalyze_comparative/dna_chunks_${input}/dna_additional_sequences*; do dnaname=$(basename $f ); rnaname="rna"${dnaname:3}; /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_for_loop.sh $dnaname $rnaname ${dnaname:3} $library $adaptor $cells $input; done

#module load R/4.1.0
bsub -q long -R "rusage[mem=16000]" -e ./$6/$4$5/log/log_mpranalyze_comparative/$8_$3_%J.e.txt -o ./$6/$4$5/log/log_mpranalyze_comparative/$8_$3_%J.o.txt Rscript /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/MPRAnalyze_for_loop_posctrl.R $1 $2 $3 $4 $5 $6 $7 $8

