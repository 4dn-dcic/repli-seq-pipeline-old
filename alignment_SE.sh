#!/bin/bash

file=$1
bowtie_genome=$2    # index name for bowtie2
nthread=$3
outprefix=$4

# log file
mapping_log_file=mapping_log.txt

# output files
bamfile=$outprefix.bam

# Mapping
bowtie2 -p $nthread -x $bowtie_genome --no-mixed --no-discordant --reorder -U $file  -S $samfile | samtools view -bS - > $bamfile 2>> $mapping_log_file
