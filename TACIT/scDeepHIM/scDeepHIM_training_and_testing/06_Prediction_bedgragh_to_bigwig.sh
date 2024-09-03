for file in *.bedGraph; do
    sorted_file="${file%.bedGraph}_sorted.bedGraph"
    bw_file="${file%.bedGraph}.bw"
    sort -k1,1 -k2,2n "$file" > "$sorted_file"
    
    bedGraphToBigWig "$sorted_file" mm10.chrom.sizes "$bw_file"
    
    echo "Processed $file to $bw_file"
done