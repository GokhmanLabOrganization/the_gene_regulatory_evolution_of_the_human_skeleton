#!/bin/bash
#
# bsub -q molgen-q -e ./$cells/$library$adaptor/log/alignBowtie_RNADNA_loop_%J.e.txt -o ./$cells/$library$adaptor/log/alignBowtie_RNADNA_loop_%J.o.txt -n 16 -R "rusage[mem=8000]" -R "span[hosts=1]" /home/labs/davidgo/Collaboration/backup/humanMPRA/scripts/alignBowtie_RNADNA_loop.sh $library $adaptor $cells $min_associations
#
threads=16
#module load bowtie2/2.3.5.1
#module load samtools/1.9 # need to specify version # Katharina 2.8.22


bowtie2-build -f ./neurons/$1$2/output/DNA_barcode_associations_$4/oligos_to_barcodes_comb_$1$2_30_$4.fasta ./$3/$1$2/output/alignments/oligos_to_barcodes_$1$2_30_$4_index # katharina 2.8.22 - build index outside of loop # katharina 5.3.23 - association output is saved in neuron folder

echo "test1"
for molecule_type in RNA DNA;
do for rep in 1 2 3;
do time bowtie2 \
    -x ./$3/$1$2/output/alignments/oligos_to_barcodes_$1$2_30_$4_index \
    -1 ./$3/$1$2/input/DNA_RNA/$1-${2}_${molecule_type}_rep${rep}_R1.fastq.gz \
    -2 ./$3/$1$2/input/DNA_RNA/$1-${2}_${molecule_type}_rep${rep}_R3.fastq.gz \
    --very-sensitive \
    --threads ${threads} \
    | samtools view \
        -@ ${threads} \
        -b \
        - \
    > ./$3/$1$2/output/alignments/Aligned_${molecule_type}_${3}_rep${rep}.bam

done;
done
