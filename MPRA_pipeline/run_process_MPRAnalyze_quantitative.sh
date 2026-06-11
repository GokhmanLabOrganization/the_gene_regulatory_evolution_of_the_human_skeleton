#!/bin/bash
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/run_process_MPRAnalyze_quantitative_%J.e.txt -o ./$cells/$library$adaptor/log/run_process_MPRAnalyze_quantitative_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_MPRAnalyze_quantitative.sh $library $adaptor $cells
# KL 20.01.23 - passing the arguments only worked when using R script insetad of R (might work with R as well but might need to give additional parameter: options "--args userargument1 userargument2") https://swcarpentry.github.io/r-novice-inflammation/05-cmdline/
module load R/4.1.0
Rscript --no-save /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/MPRAnalyze_quantitative.R $1 $2 $3
