# cells=chondrocytes
# library=L3
# adaptor=a2
# bsub -q molgen-q -R "rusage[mem=50000]" -e ./$cells/$library$adaptor/log/filter_outliers_%J.e.txt -o ./$cells/$library$adaptor/log/filter_outliers_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_filter_outliers.sh $library $adaptor $cells $min_associations

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/filter_outliers.py $1 $2 $3