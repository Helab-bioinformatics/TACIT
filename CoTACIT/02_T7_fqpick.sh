#!usr/bin/
for i in `ls *_cut27.R2.fq.gz`
do 
    base_1=$(basename $i "_cut27.R2.fq.gz")
    cd ${base_1}.R2/
    for j in `ls *.2.fq.gz`
    do base_2=$(basename $j ".2.fq.gz")
    perl ../fqpick_clhuada.pl ${base_2}.2.fq.gz ../${base_1}_cut21.R1.fq.gz ${base_2}.1.fq.gz
    done
    cd ..
done
