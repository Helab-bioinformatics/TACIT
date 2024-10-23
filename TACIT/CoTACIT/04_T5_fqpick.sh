#!usr/bin/
for i in `ls *.1.fq.gz`
do 
    base_1=$(basename $i ".1.fq.gz")
    cd ${base_1}.R1/
     for j in `ls *.1.fq.gz`
     do base_2=$(basename $j ".1.fq.gz")
     perl ../../fqpick_clhuada.pl ${base_2}.1.fq.gz ../${base_1}.2.fq.gz ${base_2}.2.fq.gz
    done
    cd ..
done
