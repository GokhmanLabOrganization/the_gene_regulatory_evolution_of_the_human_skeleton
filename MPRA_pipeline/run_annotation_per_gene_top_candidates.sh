# cells=chondrocytes
# bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./top_candidates/$cells/log/annotation_per_gene_top_candidates_%J.e.txt -o ./top_candidates/$cells/log/annotation_per_gene_top_candidates_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_annotation_per_gene_top_candidates.sh $cells

module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/annotation_per_gene_top_candidates.py $1