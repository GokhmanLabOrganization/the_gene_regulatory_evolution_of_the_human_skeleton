# bsub -q molgen-q -R "rusage[mem=1400]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/count_lines_fastq.sh

for f in *.gz; do zcat "$f" | wc -l; done