#!/bin/bash

file1=$1   # R1.fastq
file2=$2   # R2.fastq
bowtie_genome=$3    # index name for bowtie2
outprefix=$4

# log file
mapping_log_file=mapping_log.txt

# output files
bamfile=$outprefix.bam

# Mapping fastq(.gz) and making bed with values on genomic windows
# Mapping
bowtie2 -x $bowtie_genome --no-mixed --no-discordant --reorder -1 $file1 -2 $file2 -S $samfile | samtools view -bS - > $bamfile 2>> $mapping_log_file
