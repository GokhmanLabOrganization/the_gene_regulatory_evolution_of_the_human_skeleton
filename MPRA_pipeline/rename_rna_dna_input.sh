
for file in *.fastq.gz; do
    prefix="${file:0:2}"
    middle="${file:2:2}"
    suffix="${file:13:3}"
    rep="${file:16:1}"
    suffix2="${file:18:2}"
    extension=".fastq.gz"
    
    new_name="${prefix}-${middle}_${suffix}_rep${rep}_${suffix2}${extension}"
    
    mv "$file" "$new_name"
done