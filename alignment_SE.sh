#!/bin/bash

file=$1
bowtie_genome=$2    # index name for bowtie2
outprefix=$3

# log file
mapping_log_file=mapping_log.txt

# output files
bamfile=$outprefix.bam

# Mapping fastq(.gz) and making bed with values on genomic windows
# Mapping
bowtie2 -x $bowtie_genome --no-mixed --no-discordant --reorder -U $file  -S $samfile | samtools view -bS - > $bamfile 2>> $mapping_log_file
