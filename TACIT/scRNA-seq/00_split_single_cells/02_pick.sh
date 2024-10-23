for i in *.2.fq.gz
do 
base=$(basename $i ".2.fq.gz")
perl ./fqpick.pl $i  ./${i%_*}_1.fq.gz  ${base}.1.fq.gz
done

