# cells=chondrocytes
# bsub -q molgen-q -R "rusage[mem=10000]" -e ./$cells/quantitative_analysis_combined/log/get_coordinates_from_active_%J.e.txt -o ./$cells/quantitative_analysis_combined/log/get_coordinates_from_active_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_get_coordinates_from_active.sh $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/get_coordinates_from_active.py $1