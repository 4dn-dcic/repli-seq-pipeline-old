# creates an output mapping_log.txt

file=
chromsize_file=
sorted_chromsize_file=
window_size=50000
genome_windows_file=your_genome_windows.bed
bowtie_genome=/path/to/your/genome
mapping_log_file=mapping_log.txt
mapq_cut=10
samfile=${file%.fastq*}.sam
bamfile=${file%.fastq*}.bam
bedfile=${file%.fastq*}.bed
bgfile=${file%.fastq*}.bg
sort_size=5G
early_suffix=E_.bg
late_suffix=L_.bg
rt_suffix=T_.bg
merged_rt_file=merge_RT.txt


# Generating 50kb windows positions along the genome. You can change the -w value to change the windows size, and the -s value, to change the steps (to make overlapping windows for example)
sort -k1,1 -k2,2n $chromsize_file > $sorted_chromsize_file
bedtools makewindows -w $window_size -s $window_size -g $sorted_chromsize_file > $genome_windows_file

# Mapping fastq(.gz) and making bed with values on genomic windows for file in *.fastq*; do
# Mapping
(bowtie2 -x $bowtie_genome --no-mixed --no-discordant --reorder -U $file -S $samfile 2>> $mapping_log_file

# sam to bam conversion
samtools view -bSq $mapq_cut $samfile > $bamfile

# bam to bed conversion
bamToBed -i $bamfile | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S $sort_size > $bedfile
# bed line number calcul
x=`wc -l $bedfile | cut -d' ' -f 1`

# generate coverage on genomic windows
bedtools intersect -sorted -c -b $bedfile -a $genome_windows_file | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > $bgfile)

#-------

# Calculating RT
for file in *_$early_suffix; do
 paste $file ${file%$early_suffix}$late_suffix | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > ${file%$early_suffix}$rt_suffix
done

# Merging RT files
bedtools unionbedg -filler "NA" -i *$rt_suffix > $merged_rt_file
