# cells=chondrocytes
# bsub -q gsla-cpu -R "rusage[mem=10000]" -e ./$cells/quantitative_analysis_combined/log/combining_quant_results_%J.e.txt -o ./$cells/quantitative_analysis_combined/log/combining_quant_results_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_combining_quant_results.sh $cells

#module load python/3.7.9
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/combining_quant_results.py $1