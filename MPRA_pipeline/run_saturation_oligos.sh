# bsub -q molgen-q -R "rusage[mem=10000]" -e ./additional/saturation_oligos/log/sat_oligos_%J.e.txt -o ./additional/saturation_oligos/log/sat_oligos_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_saturation_oligos.sh


module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/saturation_oligos.py