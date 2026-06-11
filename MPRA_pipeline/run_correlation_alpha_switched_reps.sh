# cells=chondrocytes
# library=L1
# adaptor=a2


# bsub -q molgen-q -R "rusage[mem=10000]" -e ./additional/switching_replicates/$cells/$library$adaptor/log/correlation_alpha_%J.e.txt -o ./additional/switching_replicates/$cells/$library$adaptor/log/correlation_alpha_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_correlation_alpha_switched_reps.sh $library $adaptor $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/correlation_alpha_switched_reps.py $1 $2 $3