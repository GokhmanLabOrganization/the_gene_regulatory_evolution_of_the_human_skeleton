# bsub -q new-short -R "rusage[mem=45000]"  -e /home/labs/davidgo/Collaboration/humanMPRA/region_gene_link/chondrocytes/log/yardens_script_%J.e.txt -o /home/labs/davidgo/Collaboration/humanMPRA/region_gene_link/chondrocytes/log/yardens_script_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_yardens_script.sh

module load python/3.7.9
module load BEDTools
folder=/home/labs/davidgo/Collaboration/Lab_Tools/Map_region_to_gene
script=$folder/Map_regions_to_genes.py
config_file=/home/labs/davidgo/Collaboration/humanMPRA/region_gene_link/chondrocytes/yardens_script/config.txt
python $script $config_file