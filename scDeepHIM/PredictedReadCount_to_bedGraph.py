import numpy as np

for i in range(1, n+1):
    file_name = f"Predicted_profiles_sample_H3K4me1_{i}.txt"
    data = np.loadtxt(file_name)
    
    transposed_data = np.transpose(data)
    
    with open(f'Predicted_profiles_sample_H3K4me1_{i}.bedGraph', 'w') as fout:
        for line in transposed_data:
            chrom = str(int(line[0]))
            start = str(int(line[1]))
            end = str(int(line[2]))
            count = str(int(line[3]))
            fout.write(f'{chrom}\t{start}\t{end}\t{count}\n')
                    
for i in range(1, n+1):
    file_name = f"Predicted_profiles_H3K4me3_sample_{i}.txt"
    data = np.loadtxt(file_name)
    
    transposed_data = np.transpose(data)
    
    with open(f'Predicted_profiles_sample_H3K4me3_{i}.bedGraph', 'w') as fout:
        for line in transposed_data:
            chrom = str(int(line[0]))
            start = str(int(line[1]))
            end = str(int(line[2]))
            count = str(int(line[3]))
            fout.write(f'{chrom}\t{start}\t{end}\t{count}\n
            
for i in range(1, n+1):
    file_name = f"Predicted_profiles_H3K36me3_sample_{i}.txt"
    data = np.loadtxt(file_name)
    
    transposed_data = np.transpose(data)
    
    with open(f'Predicted_profiles_sample_H3K36me3_{i}.bedGraph', 'w') as fout:
        for line in transposed_data:
            chrom = str(int(line[0]))
            start = str(int(line[1]))
            end = str(int(line[2]))
            count = str(int(line[3]))
            fout.write(f'{chrom}\t{start}\t{end}\t{count}\n')