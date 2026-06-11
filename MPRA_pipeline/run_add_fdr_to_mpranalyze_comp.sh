# cells=chondrocytes
# bsub -q molgen-q -R "rusage[mem=10000]" -e ./$cells/comparative_analysis_combined/log/add_fdr_to_mpranalyze_comp_%J.e.txt -o ./$cells/comparative_analysis_combined/log/add_fdr_to_mpranalyze_comp_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_add_fdr_to_mpranalyze_comp.sh $cells

#module load python/3.7.9
#pip install statsmodels --user
python /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/add_fdr_to_mpranalyze_comp.py $1