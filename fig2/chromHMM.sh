###ChromHMM
wget -c http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM
# binarization
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar BinarizeBam -b 200 -f 0 -paired ./mm10.chrom.nochrM.sizes ./  cellmarkfiletable.txt binarization
# learning models
java -mx4000M -jar /media/helab/data1/min/00_reference/software/ChromHMM/ChromHMM.jar LearnModel -d 0.001 -color 0,0,255 -printposterior -p 5 -i chrhmm binarization learnmodel 10 mm10



