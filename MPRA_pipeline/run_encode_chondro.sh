# bsub -q new-short -R "rusage[mem=10000]" -e ./region_gene_link/chondrocytes/log/encode_chondro_%J.e.txt -o ./region_gene_link/chondrocytes/log/encode_chondro_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_encode_chondro.sh

module load python/3.7.9
pip install venn --user
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/encode_chondro.py