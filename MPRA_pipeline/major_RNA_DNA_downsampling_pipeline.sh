# DNA RNA - downsampling
# Created by Nadav Mishol 07-10-2024. This scripts holds all the necessary stages to perfrom a downsampling analysis for the RNA DNA reads 


# cells=chondrocytes
# library=L4
# adaptor=a1
# min_associations=2

cd /home/labs/davidgo/Collaboration/humanMPRA/


# Downsample the DNA RNA reads. NOTE: this is done non-randomly. N.M 07-10-2024 NOTE2: manually prepare the folders for each percent in the downsampling folder
for i in 1 5 25 50 70 80 90;
do bsub -q gsla-cpu -e $cells/$library$adaptor/log/RNA_DNA_downsampling/RNA_DNA_downsampling_%J.e.txt -o $cells/$library$adaptor/log/RNA_DNA_downsampling/RNA_DNA_downsampling_%J.o.txt -R "rusage[mem=140]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/splitting_RNADNA_fastqs.sh $library $adaptor $cells $i
done

for i in 1 5 25 50 70 80 90;
do bsub -q gsla-cpu -R "rusage[mem=40000]" -e ${cells}Downsampling$i/$library$adaptor/log/alignBowtie_RNADNA_loop_%J.e.txt -o ${cells}Downsampling$i/$library$adaptor/log/alignBowtie_RNADNA_loop_%J.o.txt -n 16 -R "rusage[mem=8000]" -R "span[hosts=1]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/alignBowtie_RNADNA_loop.sh $library $adaptor ${cells}Downsampling$i $min_associations
done


for i in 1 5 25 50 70 80 90;
do bsub -q gsla-cpu -R "rusage[mem=40000]" -e ${cells}Downsampling$i/$library$adaptor/log/run_process_RNA_DNA_processing_%J.e.txt -o ${cells}Downsampling$i/$library$adaptor/log/run_process_RNA_DNA_processing_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_process_RNA_DNA_processing.sh $library $adaptor ${cells}Downsampling$i $min_associations
done


for i in 1 5 25 50 70 80 90;
do bsub -q gsla-cpu -R "rusage[mem=50000]" -e ${cells}Downsampling$i/$library$adaptor/log/filter_outliers_std2_%J.e.txt -o ${cells}Downsampling$i/$library$adaptor/log/filter_outliers_std2_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_filter_outliers_std2.sh $library $adaptor ${cells}Downsampling$i $min_associations
done

#Run new and improved  correlation comparison 
bsub -q gsla-cpu -R "rusage[mem=50000]" -e additional/downsampling_corr_barcodes_oligos_Nadav/$library$adaptor/log/downsampling_statistics_%J.e.txt -o additional/downsampling_corr_barcodes_oligos_Nadav/$library$adaptor/log/downsampling_statistics_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_downsampling_statistics.sh $library $adaptor


#Run original correlation comparison 
for i in 1 5 25 50 70 80 90;
do bsub -q gsla-cpu -R "rusage[mem=50000]" -e ${cells}Downsampling$i/$library$adaptor/log/test_corr_w_downsampling_%J.e.txt -o ${cells}Downsampling$i/$library$adaptor/log/test_corr_w_downsampling_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_downsampling_corr.sh $library $adaptor ${cells}Downsampling$i $min_associations
done

# # running script for oligo coverage
# for i in 90 80 70 50 25;
# do for rep in rep1 rep2 rep3;
# do bsub -q molgen-q -R "rusage[mem=50000]" -e ./additional/rna_dna_check/$i/log/rna_dna_bcs_per_oligo_${cells}_${library}${adaptor}_${rep}_%J.e.txt -o ./additional/rna_dna_check/$i/log/rna_dna_bcs_per_oligo_${cells}_${library}${adaptor}_${rep}_%J.o.txt ./additional/rna_dna_check/run_rna_dna_bcs_per_oligo.sh $library $adaptor $cells $rep $i $min_associations;
# done;
# done

# # run script for correlation between replicates
# for i in 90 80 70 50 25;
# do bsub -q molgen-q -R "rusage[mem=50000]" -e ./additional/corr/$cells/log/rna_dna_bcs_per_oligo_${cells}_${library}${adaptor}_${i}_%J.e.txt -o ./additional/corr/$cells/log/rna_dna_bcs_per_oligo_${cells}_${library}${adaptor}_${i}_%J.o.txt /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/run_downsampling_corr.sh $library $adaptor $cells $i;
# done

# run script calculate ratio
for i in 90 80 70 50 25;
do bsub -q molgen-q -R "rusage[mem=50000]" -e ./$cells/$library$adaptor/$i/log/run_calculate_ratio_df_%J.e.txt -o ./$cells/$library$adaptor/$i/log/run_calculate_ratio_df_%J.o.txt ./run_calculate_ratio_df.sh $library $adaptor $cells $i
done

# run script quant results
for i in 90 80 70 50 25;
do bsub -q molgen-q -R "rusage[mem=140000]" -e ./$cells/$library$adaptor/$i/log/run_quant_results_%J.e.txt -o ./$cells/$library$adaptor/$i/log/run_quant_results_%J.o.txt ./run_quant_results.sh $library $adaptor $cells $i
done
