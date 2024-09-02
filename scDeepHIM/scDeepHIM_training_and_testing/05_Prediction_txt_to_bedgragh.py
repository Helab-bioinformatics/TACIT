import glob

txt_files = glob.glob("*.txt")

for txt_file in txt_files:
    output_file = txt_file.replace('.txt', '.bedGraph')
    with open(txt_file, 'r') as f:
        with open(output_file, 'w') as fout:
            for line in f:
                parts = line.strip().split('\t')
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                count = int(parts[3])
                fout.write(f'{chrom}\t{start}\t{end}\t{count}\n')
