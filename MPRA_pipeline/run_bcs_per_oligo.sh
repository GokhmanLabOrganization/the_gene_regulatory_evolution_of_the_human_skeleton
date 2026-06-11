#08-10-2024 Nadav Mishol

# bsub -q gsla-cpu  -R "rusage[mem=140000]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_bcs_per_oligo.sh $library $adaptor $cells $min_associations

# cells=chondrocytes
# library=L4
# adaptor=a1a2
# min_associations=2


python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/bcs_per_oligo.py $1 $2 $3 $4