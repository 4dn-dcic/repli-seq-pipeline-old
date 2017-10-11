#!/bin/bash

file1=$1   # R1.fastq
file2=$2   # R2.fastq
bowtie_genome=$3    # index name for bowtie2
nthread=$4
outprefix=$5

# log file
mapping_log_file=mapping_log.txt

# output files
bamfile=$outprefix.bam

# Mapping
bowtie2 -p $nthread -x $bowtie_genome --no-mixed --no-discordant --reorder -1 $file1 -2 $file2 -S $samfile | samtools view -bS - > $bamfile 2>> $mapping_log_file
