#!/bin/bash

echo "---start running ----"

ml miniconda
conda activate snakemake_hybrids

ml STAR
ml SAMtools
ml jdk
# ml UMI-tools

snakemake_job="Hy_counts"
config_file=$1
output_path=$2

bsub -q "gsla-cpu" -R "rusage[mem=50000]" -J "$snakemake_job" \
  -o "$output_path/snakemake_$(date +%Y%m%d_%H%M).o.txt" \
  -e "$output_path/snakemake_$(date +%Y%m%d_%H%M).e.txt" \
  python runSnake.py -c $config_file

echo "---end running ----"