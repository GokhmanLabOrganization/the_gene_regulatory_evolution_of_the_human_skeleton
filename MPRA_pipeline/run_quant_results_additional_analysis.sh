# cells=neurons
# library=L1
# adaptor=a1
# bsub -q molgen-q -R "rusage[mem=40000]" -e ./$cells/$library$adaptor/log/run_quant_results_additional_analysis_%J.e.txt -o ./$cells/$library$adaptor/log/run_quant_results_additional_analysis_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_quant_results_additional_analysis.sh $library $adaptor $cells


module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/quant_results_additional_analysis.py $1 $2 $3