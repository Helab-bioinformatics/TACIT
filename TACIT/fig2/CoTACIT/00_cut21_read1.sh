for i in `ls *_1.fq.gz`; do base=$(basename $i "_1.fq.gz"); cutadapt  -u 21 -o 01_cut/${base}_cut21.R1.fq.gz ${base}_1.fq.gz; done
