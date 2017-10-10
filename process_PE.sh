file1=$1   # R1.fastq
file2=$2   # R2.fastq
bowtie_genome=$3    # index name for bowtie2
chromsize_file=$4
outprefix=$5

# parameters
window_size=50000
mapq_cut=10
sort_size=5G

# intermediate files, data-independent
sorted_chromsize_file=$chromsize_file.sorted
genome_windows_file=$chromsize_file.windows.bed

# log file
mapping_log_file=mapping_log.txt

# output files
samfile=$outprefix.sam
bamfile=$outprefix.bam
bedfile=$outprefix.bed
bgfile=$outprefix.bg

# Generating the 50kb windows positions along the genome. You can change the -w value to change the windows size, and the -s value, to change the steps (to make overlapping windows for example)
sort -k1,1 -k2,2n $chromsize_file > $sorted_chromsize_file
bedtools makewindows -w $window_size -s $window_size -g $sorted_chromsize_file > $genome_windows_file


# Mapping fastq(.gz) and making bed with values on genomic windows
# Mapping
(bowtie2 -x $bowtie_genome --no-mixed --no-discordant --reorder -1 $file1 -2 $file2 -S $samfile 2>> $mapping_log_file
# sam to bam conversion
samtools view -bSq $mapq_cut $samfile > $bamfile
# bam to bed conversion
bamToBed -i $bamfile | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S $sort_size > $bedfile
# bed line number calcul
x=`wc -l $bedfile | cut -d' ' -f 1`
# generate coverage on genomic windows
bedtools intersect -sorted -c -b $bedfile -a $genome_windows_file | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > $bgfile)

