#!/bin/bash

input_bamfile=$1   # unfiltered bam file
chromsize_file=$2
outprefix=$3

# parameters
window_size=50000
mapq_cut=10
sort_size=5G

# intermediate files, data-independent
sorted_chromsize_file=$chromsize_file.sorted
genome_windows_file=$chromsize_file.windows.bed

# output files
bamfile=$outprefix.bam  # filtered bam file
bedfile=$outprefix.bed
bgfile=$outprefix.bg

# Generating the 50kb windows positions along the genome. You can change the -w value to change the windows size, and the -s value, to change the steps (to make overlapping windows for example)
sort -k1,1 -k2,2n $chromsize_file > $sorted_chromsize_file
bedtools makewindows -w $window_size -s $window_size -g $sorted_chromsize_file > $genome_windows_file

# filtering
samtools view -bq $mapq_cut $input_bamfile > $bamfile
# bam to bed conversion
bamToBed -i $bamfile | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S $sort_size > $bedfile
# bed line number calcul
x=`wc -l $bedfile | cut -d' ' -f 1`
# generate coverage on genomic windows
bedtools intersect -sorted -c -b $bedfile -a $genome_windows_file | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > $bgfile)
