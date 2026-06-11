# cells=chondrocytes

# bsub -q gsla-cpu -R "rusage[mem=50000]" -e ./additional/analyze_controls/$cells/log/analyze_controls_%J.e.txt -o ./additional/analyze_controls/$cells/log/analyze_controls_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_analyze_controls.sh $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/analyze_controls.py $1