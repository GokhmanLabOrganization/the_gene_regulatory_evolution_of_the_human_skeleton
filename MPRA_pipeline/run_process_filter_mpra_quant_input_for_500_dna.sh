# cells=neurons
# library=L1
# adaptor=a2
# bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/filter_mpra_quant_input_for_500_dna_%J.e.txt -o ./$cells/$library$adaptor/log/filter_mpra_quant_input_for_500_dna_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_filter_mpra_quant_input_for_500_dna.sh $library $adaptor $cells


module load python/3.7.9

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/filter_mpra_quant_input_for_500_dna.py $1 $2 $3