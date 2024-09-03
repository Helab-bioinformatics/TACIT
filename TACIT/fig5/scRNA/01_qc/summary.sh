for i in `ls ./00_raw_data/*.1.fq.gz`
do 
base=$(basename $i ".1.fq.gz")
zcat ./00_raw_data/${base}.1.fq.gz > ./00_raw_data/${base}.1.fq
wc -l ./00_raw_data/${base}.1.fq >> QC.xls
done

rm ./00_raw_data/*.1.fq

Rscript summary.R