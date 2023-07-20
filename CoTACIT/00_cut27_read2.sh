for i in `ls *_2.fq.gz`; do base=$(basename $i "_2.fq.gz"); cutadapt  -u 27 -o 01_cut/${base}_cut27.R2.fq.gz ${base}_2.fq.gz; done
