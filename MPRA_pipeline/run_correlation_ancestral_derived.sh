# cells=neurons
# library=L2
# adaptor=a2
# bsub -q new-short -R "rusage[mem=10000]" -e ./$cells/$library$adaptor/log/correlation_ancestral_derived_%J.e.txt -o ./$cells/$library$adaptor/log/correlation_ancestral_derived_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_correlation_ancestral_derived.sh $library $adaptor $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/correlation_ancestral_derived.py $1 $2 $3