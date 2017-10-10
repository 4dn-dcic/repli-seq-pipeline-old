#!/bin/bash
## Usage: merge_rt.sh outdir in1.bg in2.bg in3.bg ...

outdir=$1
shift

INFILESTR=${@}  # input files as a string
INFILES=(${INFILESTR// / })  # input files as an array

rt_files=$INFILESTR  # array of files, as a space-delimited string
merged_rt_file=$outdir/merged_rt.bg

# Merging RT files
bedtools unionbedg -filler "NA" -i $rt_files > $merged_rt_file
