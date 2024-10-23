for i in `ls *.1.fq.gz`; do base=$(basename $i ".1.fq.gz"); mkdir ${base}.R1/; perl ../splitbarcode0bp.pl ${base}.1.fq.gz 0 AAAGAA CAAGTC GAAGGG ${base}.R1/T5_${base}_1.1.fq.gz ${base}.R1/T5_${base}_2.1.fq.gz ${base}.R1/T5_${base}_3.1.fq.gz ${base}.R1/T5_${base}_noid.1.fq.gz; done



