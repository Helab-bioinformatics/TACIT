# Data preprocessing, initiated from bam files
# Step1: convert bam files to bigwig files
for i in *_rmdup_picard.bam
do
samtools index $i
done

for i in ./*_rmdup_picard.bam
do
base=$(basename $i "_rmdup_picard.bam")
num1=10000000
num2="$(samtools view -c  $i  2>&1 )"
res=$(printf "%.5f" `echo "scale=5;$num1/$num2"|bc`)
bamCoverage --scaleFactor  $res -b  $i  -o  ./${base}.bw
done

# Step2: compute 10kb signals
multiBigwigSummary BED-file -b *bw  -out mm10_genome10kb_signals.results.npz --BED   mm10.10K.windows.bed --outRawCounts mm10_genome10kb_signals_readCount.tab 
