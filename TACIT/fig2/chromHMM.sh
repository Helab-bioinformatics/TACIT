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

