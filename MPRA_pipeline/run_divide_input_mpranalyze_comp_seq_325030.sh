# bsub -q new-short -R "rusage[mem=40000]" -e ./chondrocytes/L3a2/log/divide_input_mpranalyze_comp_seq_325030_%J.e.txt -o ./chondrocytes/L3a2/log/divide_input_mpranalyze_comp_seq_325030_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_divide_input_mpranalyze_comp_seq_325030.sh


module load python/3.7.9

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/divide_input_mpranalyze_comp_seq_325030.py