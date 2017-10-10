# Generating 50kb windows positions along the genome. You can change the -w value to change the windows size, and the -s value, to change the steps (to make overlapping windows for example)
sort -k1,1 -k2,2n your_genome.chrom.sizes > your_genome_sorted.chrom.sizes
bedtools makewindows -w 50000 -s 50000 -g your_genome_sorted.chrom.sizes > your_genome_windows.bed

# Mapping fastq(.gz) and making bed with values on genomic windows for file in *.fastq*; do
# Mapping
(bowtie2 -x /path/to/your/genome --no-mixed --no-discordant --reorder -U $file -S ${file%.fastq*}.sam 2>> mapping_log.txt

# sam to bam conversion
samtools view -bSq 10 ${file%.fastq*}.sam > ${file%.fastq*}.bam

# bam to bed conversion
bamToBed -i ${file%.fastq*}.bam | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S 5G > ${file%.fastq*}.bed
# bed line number calcul
x=`wc -l ${file%.fastq*}.bed | cut -d' ' -f 1`

# generate coverage on genomic windows
bedtools intersect -sorted -c -b ${file%.fastq*}.bed -a your_genome_windows.bed | awk -vx=$x
'{print $1,$2,$3,$4*1e+06/x}' OFS='\t' > ${file%.fastq*}.bg) &
done
wait

# Calculating RT
for file in *_E_.bg; do
 paste $file ${file%E_.bg}L_.bg | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > ${file%E_.bg}T_.bg
done

# Merging RT files
bedtools unionbedg -filler "NA" -i *T_.bg > merge_RT.txt
