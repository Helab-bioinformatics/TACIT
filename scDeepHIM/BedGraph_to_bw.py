import os
for i in range(1, n+1):
    os.system(f'sort -k1,1 -k2,2n Predicted_profiles_sample_H3K36me3_{i}.bedGraph > Predicted_profiles_sample_H3K36me3_{i}_sorted.bedGraph')
    os.system(f'sort -k1,1 -k2,2n Predicted_profiles_sample_H3K4me3_{i}.bedGraph > Predicted_profiles_sample_H3K4me3_{i}_sorted.bedGraph')
    os.system(f'sort -k1,1 -k2,2n Predicted_profiles_sample_H3K4me1_{i}.bedGraph > Predicted_profiles_sample_H3K4me1_{i}_sorted.bedGraph')
    os.system(f'bedGraphToBigWig Predicted_profiles_sample_H3K36me3_{i}_sorted.bedGraph mm10.chrom.sizes Predicted_profiles_sample_H3K36me3_{i}.bw')
    os.system(f'bedGraphToBigWig Predicted_profiles_sample_H3K4me3_{i}_sorted.bedGraph mm10.chrom.sizes Predicted_profiles_sample_H3K4me3_{i}.bw')
    os.system(f'bedGraphToBigWig Predicted_profiles_sample_H3K4me1_{i}_sorted.bedGraph mm10.chrom.sizes Predicted_profiles_sample_H3K4me1_{i}.bw')
