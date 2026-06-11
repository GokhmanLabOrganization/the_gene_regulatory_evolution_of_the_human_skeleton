# cells=chondrocytes
# bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./$cells/comparative_analysis_combined/log/get_coordinates_L1L2L3L4_%J.e.txt -o ./$cells/comparative_analysis_combined/log/get_coordinates_L1L2L3L4a1_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_get_coordinates_L1L2L3L4a1.sh $cells

module load python/3.7.9
python3 /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/get_coordinates_L1L2L3L4a1.py $1