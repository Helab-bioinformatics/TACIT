for i in `ls *_1.fq.gz`; do base=$(basename $i "_1.fq.gz"); nohup bowtie2 -x /media/helab/data1/00_public/database/genomes/mm10_bowtie2/mm10 -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz --dovetail --very-sensitive-local --no-unal --no-mixed --no-discordant -X 2000 -S ${base}.mm10.sam 2> ${base}.bowtie2.align.log -p10; done

for i in `ls *.sam`; do base=$(basename $i ".sam"); samtools view -hbS -q 30 $i | samtools sort -T ${base}.mm10 - > ${base}_unique.bam; done
for i in `ls *.bowtie2.align.log`; do base=$(basename $i ".bowtie2.align.log");cat ${base}.bowtie2.align.log | awk '/overall/{print $1}' >>Mapping.txt; echo ${base} >>MappingID.txt; done

for i in `ls *_unique.bam`; do base=$(basename $i "_unique.bam");java -jar /media/helab/data1/00_public/software/picard-tools-2.2.4/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$i o=${base}_rmdup.bam M=${base}_picard.txt;done
for i in `ls *_unique.bam`; do base=$(basename $i "_unique.bam"); bamtools stats -in ${base}_unique.bam > ${base}.sam.uni.totalreads.log; done

for i in `ls *.sam.uni.totalreads.log`; do base=$(basename $i ".sam.uni.totalreads.log");cat ${base}.sam.uni.totalreads.log | awk '/Total/{print $3}' >>Unimap.txt; echo ${base}>>UnimapID.txt; done
for i in `ls *_rmdup.bam`; do base=$(basename $i "_rmdup.bam"); bamtools stats -in ${base}_rmdup.bam > ${base}.sam.rmdup.totalreads.log; done
for i in `ls *.sam.rmdup.totalreads.log`; do base=$(basename $i ".sam.rmdup.totalreads.log"); cat ${base}.sam.rmdup.totalreads.log | awk '/Total/{print $3}' >>rmdup.txt; echo ${base} >>rmdupID.txt; done

