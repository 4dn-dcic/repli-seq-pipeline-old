rt_files=  # array of files
merged_rt_file=

# Merging RT files
bedtools unionbedg -filler "NA" -i $rt_files > $merged_rt_file
