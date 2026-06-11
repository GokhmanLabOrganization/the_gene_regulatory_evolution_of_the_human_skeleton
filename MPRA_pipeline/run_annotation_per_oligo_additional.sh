# cells=chondrocytes
# bsub -q molgen-q -R "rusage[mem=10000]" -e ./top_candidates/$cells/log/annotation_per_oligo_additional_%J.e.txt -o ./top_candidates/$cells/log/annotation_per_oligo_additional_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_annotation_per_oligo_additional.sh $cells

module load python/3.7.0
module load BEDTools
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/annotation_per_oligo_additional.py $1