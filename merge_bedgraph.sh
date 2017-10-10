## (Optional) Merge your bedgraph files for further analysis. 
cd path/to/files/ 
#  Merge the loess smoothed bedgraph files:
bedtools unionbedg -filler "NA" -i *Loess.bedGraph > merge_Loess_norm_RT.txt 
