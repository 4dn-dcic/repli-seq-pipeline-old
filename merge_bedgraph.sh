#!/bin/bash
## (Optional) Merge your bedgraph files for further analysis. 
## usage: merge_bedgraph.sh outdir in1.bg in2.bg in3.bg ...
## output file merge_Loess_norm_RT.txt will be created in outdir.

outdir=$1
shift

INFILESTR=${@}  # input files as a string
INFILES=(${INFILESTR// / })  # input files as an array

# input files, loess smoothed bedgraph files
loess_bg_files=$INFILESTR  # array of files, as a space-delimited string.

# output file name
outfile=$outdir/merge_Loess_norm_RT.txt 

#  Merge the loess smoothed bedgraph files:
bedtools unionbedg -filler "NA" -i $loess_bg_files > $outfile 
