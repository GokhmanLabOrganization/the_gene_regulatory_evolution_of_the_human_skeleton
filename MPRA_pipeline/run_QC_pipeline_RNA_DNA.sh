

# bsub -q gsla-cpu -R "rusage[mem=50000]" -e ./additional/QC_pipeline/log/QC_pipeline_%J.e.txt -o ./additional/QC_pipeline/log/QC_pipeline_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_QC_pipeline_RNA_DNA.sh

python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/QC_pipeline_RNA_DNA.py $1 $2 $3