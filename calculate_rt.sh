#!/bin/bash

# input files
early_file=$1
late_file=$2

# output file
rt_file=$3

# Calculating RT
paste $early_file $late_file | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > $rt_file
