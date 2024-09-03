###ChromHMM
wget -c http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM
#第一步把所有bam文件binned
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar BinarizeBam -b 200 -f 0 -paired /media/helab/data1/min/00_reference/mm10.chrom.sizes ./  cellmarkfiletable.txt binarization
#第二步统计binned文件，得到染色质富集模式
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar LearnModel -d 0.001 -color 0,0,255 -printposterior -p 5 -i chrhmm binarization learnmodel 10 mm10

java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar BinarizeBam -b 200 -f 0 -paired  ../mm10.chrom.nochrM.sizes ./  cellmarkfiletable.txt binarization_nochrM

#### makeSegamentaion
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar MakeSegmentation  -printposterior -i chrhmm /media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/ICM/learnmodel_nochrM/model_10_chrhmm.txt binarization_nochrM learnmodel_ICM 


##### BinarizeSignal模式 ######
# 参见/media/helab/data3/Moue_public/20240426_integration/paired/2cell
# 需要先准备好signal文件
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar BinarizeSignal signaldir outputdir
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar LearnModel -d 0.001 -color 0,0,255 -printposterior -p 5 -i chrhmm outputdir learnmodel 10 mm10


##### 或者直接binarization，自己准备好相关的binarization ######
# 参见/media/helab/data3/Moue_public/20240426_integration/paired/2cell
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar LearnModel -d 0.001 -color 0,0,255 -printposterior -p 5 -i chrhmm binarization learnmodel 10 mm10


################################# 得到指定区域的chromatin stats的信息 #################################
cd 01_ICM_TE_enhancer/
cat ICM_mm10_chrhmm_enhancer_strong.bed ICM_mm10_chrhmm_enhancer_weak.bed >ICM_mm10_chrhmm_enhancer.bed
cat TE_mm10_chrhmm_enhancer_strong.bed TE_mm10_chrhmm_enhancer_weak.bed >TE_mm10_chrhmm_enhancer.bed
mergePeaks -prefix mytest TE_mm10_chrhmm_enhancer.bed ICM_mm10_chrhmm_enhancer.bed

sed -i '1d' ICM_only_active.bed
awk '{print $2, $3, $4}' ICM_only_active.bed > ICM_only_active.2.bed
awk 'BEGIN{ FS=" ";OFS="\t" }{ print $1,$2,$3}' ICM_only_active.2.bed > ICM_only_active.3.bed

for i in `ls *bed`;do base=$(basename $i ".bed"); sed -i '1d' $i;awk '{print $2, $3, $4}' $i > ${base}.2.bed;awk 'BEGIN{ FS=" ";OFS="\t" }{ print $1,$2,$3}' ${base}.2.bed >${base}.3.bed;done


####把输出的bed分bin，用于不同组之间比较
cp /media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/00_synthetic_cells/05_states_concatenated/learnmodel_nochrM_10/*segments.bed ./

for i in ./*.bed
do
base=$(basename $i ".bed")
bedtools makewindows -b $i -w 200 -i src > ./${base}_200.txt
done



##删除行末制表符
awk '{gsub(/^[\t ]*|[\t ]*$/,"");print}' ICM_only_enhancer_200.txt >ICM_only_enhancer_200.2.txt

for i in ./*_segments_200.txt
do
base=$(basename $i "_segments_200.txt")
bedtools intersect -a $i -b TE_only_active_200.txt -wa >${base}_inte_TE_only_active.txt
done


for i in ./*.bed
do
base=$(basename $i ".bed")
bedtools makewindows -b $i -w 200 -i src > ./${base}_200.txt
cat ${base}_200.txt | wc -l >>line.txt
echo ${base} >>ID.txt
done