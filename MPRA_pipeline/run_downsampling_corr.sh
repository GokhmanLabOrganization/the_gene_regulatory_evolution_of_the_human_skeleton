#Jobs submission:
#bsub -q gsla-cpu -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/log/test_corr_w_downsampling_%J.e.txt -o ./$cells/$library$adaptor/log/test_corr_w_downsampling_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_downsampling_corr.sh $library $adaptor $cells

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/downsampling_corr.py $1 $2 $3 $4