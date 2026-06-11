# cells=neurons
# library=L3
# adaptor=a1
# bsub -q new-short -R "rusage[mem=50000]" -e ./$cells/$library$adaptor/log/filter_outliers_std2_%J.e.txt -o ./$cells/$library$adaptor/log/filter_outliers_std2_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_filter_outliers_std2.sh $library $adaptor $cells $min_associations

#module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/filter_outliers_std2.py $1 $2 $3