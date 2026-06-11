# bsub -q molgen-q -R "rusage[mem=10000]" -e ./additional/chengyu/log/compare_neurons_to_chengyu_%J.e.txt -o ./additional/chengyu/log/compare_neurons_to_chengyu_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_compare_neurons_to_chengyu.sh

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/compare_neurons_to_chengyu.py