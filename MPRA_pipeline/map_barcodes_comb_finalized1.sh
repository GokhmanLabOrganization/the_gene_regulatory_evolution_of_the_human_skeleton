#!/bin/bash
threads=16
module load bowtie2/2.3.5.1
module load samtools/1.9 # need to specify version # Katharina 27.7.22

# Building an index for the reference sequence
# bowtie2-build -f ./oligo_fasta/$1$2.fasta ./oligo_fasta/$1$2_index
# -x is the index name
# -1 and -2 are the paired seq files
# very-sensitive is a collection of default parameters
# samtools view converts the output to a bam file

time bowtie2 \
    -x ./oligo_fasta/$1$2_index \
    -1 ./$3/$1$2/input/DNA_barcode/$1-$2_R1.fastq.gz \
    -2 ./$3/$1$2/input/DNA_barcode/$1-$2_R3.fastq.gz \
    --very-sensitive \
    --threads ${threads} \
    | samtools view \
        -@ ${threads} \
        -b \
        - \
    > ./$3/$1$2/output/Aligned_DNA_barcodes_comb.bam

# Suggestion: create a sorted BAM file and an indexed BAM file
#samtools sort -m 4G  -o Aligned_DNA_barcodes_comb_sorted.bam Aligned_DNA_barcodes_comb.bam
#samtools index Aligned_DNA_barcodes_comb_sorted.bam
