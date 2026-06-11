# cells=chondrocytes
# bsub -q molgen-q -R "rusage[mem=10000]" -e ./$cells/comparative_analysis_combined/log/combine_quant_comp_%J.e.txt -o ./$cells/comparative_analysis_combined/log/combine_quant_comp_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_combine_quant_comp.sh $cells

module load python/3.7.9
python3 /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/combine_quant_comp.py $1