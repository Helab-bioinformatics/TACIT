###download chromHMM
wget -c http://compbio.mit.edu/ChromHMM/ChromHMM.zip
unzip ChromHMM

### binarization
java -mx4000M -jar ./ChromHMM/ChromHMM.jar BinarizeBam -b 200 -f 0 -paired ./mm10.chrom.nochrM.sizes ./  cellmarkfiletable.txt binarization

### learnmodel for each synthetic cell
java -mx4000M -jar ./ChromHMM/ChromHMM.jar LearnModel -d 0.001 -color 0,0,255 -printposterior -p 15 -i chrhmm binarization learnmodel 12 mm10

### makeSegmentation with ICM_TE_model
java -mx4000M -jar ./ChromHMM/ChromHMM.jar MakeSegmentation  -printposterior -i ./ICM_TE_model_12_chrhmm.txt binarization learnmodel

