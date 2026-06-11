
# library=L1
# adaptor=a3
# bsub -q molgen-q -R "rusage[mem=40000]" -e ./additional/chondro_neuron_comparison_thesis/log/run_chondro_neuron_compare_%J.e.txt -o ./additional/chondro_neuron_comparison_thesis/log/run_chondro_neuron_compare_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_chondro_neuron_compare_thesis.sh $library $adaptor


module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/chondro_neuron_compare_thesis.py $1 $2