early_file=
late_file=
rt_file=

# Calculating RT
paste $early_file $late_file | awk '{if($8 != 0 && $4 != 0){print $1,$2,$3,log($4/$8)/log(2)}}' OFS='\t' > $rt_file
