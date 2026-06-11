#08-10-2024 Nadav Mishol
#This script downsamplese RNA DNA data. It does so by treating downsampled data as new cell types, so downstream analysis could use the original major pipeline. 

# bsub -q gsla-cpu -e $cells/$library$adaptor/log/RNA_DNA_downsampling/alignBowtie_RNADNA_loop_%J.e.txt -o $cells/$library$adaptor/log/RNA_DNA_downsampling/alignBowtie_RNADNA_loop_%J.o.txt -R "rusage[mem=140000]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/splitting_RNADNA_fastqs.sh $library $adaptor $cells


# cells=chondrocytes
# library=L4
# adaptor=a1a2
# min_associations=2


library=$1
adaptor=$2
cells=$3
sample=$4


cd /home/labs/davidgo/Collaboration/humanMPRA/

for rep in rep1 rep2 rep3;
do for mol in DNA RNA;
do total_lines=$(zcat /home/labs/davidgo/Collaboration/humanMPRA/$cells/$library$adaptor/input/DNA_RNA/$library-${adaptor}_${mol}_${rep}_R1.fastq.gz | wc -l);
echo $total_lines;
for i in sample;
do ans=$(($(($(($((total_lines / 4)) * $sample ))/100))*4));
echo $ans;
for j in 1 2 3;
do zcat $cells/$library$adaptor/input/DNA_RNA/$library-${adaptor}_${mol}_${rep}_R${j}.fastq.gz | head -n $ans | gzip > ${cells}Downsampling$sample/$library$adaptor/input/DNA_RNA/$library-${adaptor}_${mol}_${rep}_R${j}.fastq.gz;
echo $sample $j $rep $mol;
done;
done;
done;
done
