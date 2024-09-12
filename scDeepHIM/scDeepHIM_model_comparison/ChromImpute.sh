#!/bin/sh

workpath=/ChromImpute/
ChromImpute=/ChromImpute/ChromImpute.jar
bedgraph=/ChromImpute/INPUTDATADIR/bedgraph
inputinfofile_all=/ChromImpute/inputinfofile_bedgraph.txt
inputinfofile=/ChromImpute/inputinfofile_bedgraph_train.txt
chromsize=/ChromImpute/mm10.chrom.sizes
CONVERTEDDATADIR=${workpath}/CONVERTEDDATADIR
DISTANCEDIR=${workpath}/DISTANCEDIR
PREDICTORDIR=${workpath}/PREDICTORDIR
TRAINDATA=${workpath}/TRAINDATA
OUTPUTDATA=${workpath}/OUTPUTDATA

java -mx100000M -jar ${ChromImpute} Convert -r 10000 ${bedgraph} ${inputinfofile_all} $chromsize $CONVERTEDDATADIR
java -mx100000M -jar ${ChromImpute} ComputeGlobalDist -r 10000 $CONVERTEDDATADIR $inputinfofile $chromsize $DISTANCEDIR

for mark in H3K36me3 H3K4me1 H3K4me3
do
echo $mark
java -mx50000M -jar ${ChromImpute} GenerateTrainData -r 10000 $CONVERTEDDATADIR $DISTANCEDIR $inputinfofile $chromsize $TRAINDATA ${mark}

for test_sample in 2cell 4cell 8cell Blastocyte Morula 
do
echo ${test_sample}
java -mx50000M -jar ${ChromImpute} Train $TRAINDATA $inputinfofile $PREDICTORDIR $test_sample $mark
java -mx50000M -jar ${ChromImpute} Apply -r 10000 $CONVERTEDDATADIR $DISTANCEDIR $PREDICTORDIR $inputinfofile $chromsize $OUTPUTDATA $test_sample $mark
done
done
