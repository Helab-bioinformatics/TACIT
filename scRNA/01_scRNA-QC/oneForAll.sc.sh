mkdir 00_raw_data
mkdir 01_cutadapt
mkdir 02_mapping
mkdir 03_count
mkdir 04_umi/

rm *_noid.1.fq.gz
rm *_noid.2.fq.gz
mv *.1.fq.gz ./00_raw_data
mv *.2.fq.gz ./00_raw_data

#### add UMI to the R1
for i in `ls ./00_raw_data/*.1.fq.gz`
do 
base=$(basename $i ".1.fq.gz")
perl add_UMI.pl ./00_raw_data/${base}.1.fq.gz ./00_raw_data/${base}.2.fq.gz ./00_raw_data/${base}.umi.fq.gz
done

#### cutadapt
for i in `ls ./00_raw_data/*.umi.fq.gz`
do 
base=$(basename $i ".umi.fq.gz")
cutadapt -a CTGTCTCTTATACACATCTGACGCTGCCGACGA -q 20 -O 10  --trim-n  -m 30  --max-n 0.1  -o  ./01_cutadapt/${base}.umi.trim.fq.gz   ./00_raw_data/${base}.umi.fq.gz
done

#### mapping by Hisat2 
#set cores
for i in `ls ./01_cutadapt/*.umi.trim.fq.gz`
do 
base=$(basename $i ".umi.trim.fq.gz")
hisat2 -x /media/helab/data1/00_public/database/Hisat2/hisat2_ref/mm10/mm10 -p 10 -U ./01_cutadapt/${base}.umi.trim.fq.gz -S ./02_mapping/${base}_align.mm10.sam 2>./02_mapping/${base}.mm10.align.log
#hisat2 -x /media/helab/data1/00_public/database/Hisat2/hisat2_ref/hg19/hg19 -p 20 -U ./01_cutadapt/${base}.umi.trim.fq.gz -S ./02_mapping/${base}_align.hg19.sam 2>./02_mapping/${base}.hg19.align.log
done

#### HTSeq count
for i in `ls ./02_mapping/*_align.mm10.sam`
do 
base=$(basename $i "_align.mm10.sam")
python2 /media/helab/data1/00_public/software/HTSeq/count.transcript.py  -s no -f sam -a 10 -o ./03_count/${base}_align.mm10.umi.sort.sam ./02_mapping/${base}_align.mm10.sam /media/helab/data1/00_public/database/Hisat2/hisat2_gtf/mm10_Refseq.transcript.uniq2.gtf >./03_count/${base}.dexseq.txt  
done

####  stat-umi
for i in `ls ./03_count/*_align.mm10.umi.sort.sam`
do 
base=$(basename $i "_align.mm10.umi.sort.sam")
python2 ./02_UMI_HTseq.py ./03_count/${base}_align.mm10.umi.sort.sam  ./04_umi/${base}.umi.txt
done


#### calculate counts and TPM of each cell #¼ÇµÃÐÞ¸Ämm10_tpm.pyÀïµÄÂ·¾¶
for i in `ls ./03_count/*_align.mm10.umi.sort.sam`
do 
base=$(basename $i "_align.mm10.umi.sort.sam")
python2 mm10_tpm.py ./04_umi/${base}.umi.txt ./04_umi/${base}.umi_tpm_gene.xls
done

#### generate count and TPM matrix
cd ./04_umi
Rscript ../merge.HTseq_count.R umi_tpm_gene.xls  geneid ./ HTSeq.count.txt &
Rscript ../merge.HTseq_TPM.R umi_tpm_gene.xls  geneid ./ HTSeq.TPM.txt &

#### summary
cd ..
sh summary.sh
