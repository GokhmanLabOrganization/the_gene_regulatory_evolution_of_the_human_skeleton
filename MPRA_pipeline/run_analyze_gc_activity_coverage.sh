# cells=neurons
# library=L4
# adaptor=a1

# bsub -q gsla-cpu -R "rusage[mem=50000]" -e ./additional/GC_activity_coverage/$cells/log/analyze_gc_activity_coverage_%J.e.txt -o ./additional/GC_activity_coverage/$cells/log/analyze_gc_activity_coverage_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_analyze_gc_activity_coverage.sh $library $adaptor $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/analyze_gc_activity_coverage.py $1 $2 $3