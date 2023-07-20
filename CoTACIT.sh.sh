

vim 01_T7_split.sh
vim 03_T5_split.sh
#对于华大平台，需要把所有fqpick.pl替换成fqpick_clhuada.pl
vim oneForAll.spilt.sh
vim 02_T7_fqpick.sh
vim 04_T5_fqpick.sh


find ./ -name "*raw_2_2.*" -exec cp {} ./h3k9me3 \;
for i in `ls *_1_2.1.fq.gz`; do base=$(basename $i "_1_2.1.fq.gz"); cat $i ${base}_2_1.1.fq.gz >${base}_k4k27ac.1.fq.gz; cat ${base}_1_2.2.fq.gz ${base}_2_1.2.fq.gz >${base}_k4k27ac.2.fq.gz;done


