# repli-seq-pipeline

* 1. `alignment_SE.sh` or `alignment_PE.sh` : takes fastq files, produces bam file
* 2. `post_alignment.sh` : takes bam file, produces bg file
* 3. `calculate_rt.sh` : takes early and late bg files, produces rt bg file
* 4. `merge_rt.sh` : takes in multiple rg bg files, produces merged rt bg file
* 5. `post_process.R` : takes merged rt bg file, produces normalized rt bg file
* 6. (optionally) `merge_bedgraph.sh` : takes multiple normalized rt bg files, produces merged normalized rt bg file

### Requirement
* bowtie2
* samtools
* bedtools
* R
  * library: preprocessCore
  
