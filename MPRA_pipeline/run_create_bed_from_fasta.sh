# cells=neurons
# library=L1
# adaptor=a1


# bsub -q molgen-q -R "rusage[mem=10000]" -e ./$cells/$library$adaptor/log/run_create_bed_from_fasta_%J.e.txt -o ./$cells/$library$adaptor/log/run_create_bed_from_fasta_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_create_bed_from_fasta.sh $library $adaptor $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/create_bed_from_fasta.py $1 $2 $3