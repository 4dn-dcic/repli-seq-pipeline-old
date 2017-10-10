#!/bin/bash
## (Optional) Merge your bedgraph files for further analysis. 

# input files, loess smoothed bedgraph files
loess_bg_files=@_  # array of files

# output file name
outfile=merge_Loess_norm_RT.txt 

cd path/to/files/ 
#  Merge the loess smoothed bedgraph files:
bedtools unionbedg -filler "NA" -i $loess_bg_files > $outfile 
