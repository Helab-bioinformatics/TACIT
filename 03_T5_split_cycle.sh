#!usr/bin/
for i in `ls *_cut27.R2.fq.gz`
do
    base_1=$(basename $i "_cut27.R2.fq.gz")
    cd ${base_1}.R2/
    sh ../03_T5_split.sh
    cd ..
done


