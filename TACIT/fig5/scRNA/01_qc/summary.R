library(gtools)
qc <- read.csv("QC.xls",header = F,sep = " ")
qc <- qc[,c(2,1)]
qc[,2] <- qc[,2]/4
rownames(qc) <- as.character(qc[,1])
rownames(qc) <- gsub("./00_raw_data/","",rownames(qc))
rownames(qc) <- gsub(".1.fq","",rownames(qc))
qc[,1] <- rownames(qc)
qc <- qc[mixedsort(rownames(qc),decreasing=F),]#sort
colnames(qc) <- c("sampleName","rawFrag")


#path1 <-  "./00_raw_datafastqc.results"
#path2 <-  "./01_cutadaptfastqc.cutadapt.results"
path3 <- "./02_mapping"
#path4 <- "./"
#path5 <- "./"

#df1 <- read.table(paste0(path1,"summary.fastqc.sort.txt"))
#df1 <- df1[seq(from=1,to=nrow(df1),by=2),]
#rownames(df1) <- as.character(df1[,1])
#df1 <- df1[mixedsort(as.character(df1[,1])),]


#df2 <- read.table(paste0(path2,"summary.fastqc.cutadapt.sort.txt"))
#df2 <- df2[seq(from=1,to=nrow(df2),by=2),]
#rownames(df2) <- as.character(df2[,1])
#df2 <- df2[mixedsort(as.character(df2[,1])),]

library(gtools)
read.files_path.align <- dir(path3, pattern = "align.log$" ,full.names = T)
read.files_path.align <- mixedsort(read.files_path.align,decreasing=F)
df3 <- c()
for (file in read.files_path.align){
  data3 <- read.delim(file)
  #df3 <- c(df3,substr(as.character(data3[14,1]),1,6))
  df3 <- c(df3,substr(as.character(data3[5,1]),1,6))
}
read.files_path.align[1:5]
#write.table(df3,"df3.txt")
qc$mapping_rate <- df3


count <- read.table("./04_umi/HTSeq.count.txt",header=T,row.names=1,sep="\t")
#count <- count[,-ncol(count)]
#sampleName <- read.table("sampleName.txt",header=F,sep="\t")
#colnames(count)<- sampleName$V1
#delete.index <- c(grep("feature",rownames(count)),grep("ambiguous",rownames(count)),grep("alignment",rownames(count)),grep("aQual",rownames(count)),grep("aligned",rownames(count))) #delete __no_feature	__ambiguous	__too_low_aQual	__not_aligned __alignment_not_unique
#count <- count[-delete.index,]
colnames(count)[1:5]
####R读文件名的时候，有时会把"-"读成".",会影响文件名的sort，因此要检查文件名
colnames(count) <- gsub("[.]","-",colnames(count))
count <- count[,mixedsort(colnames(count),decreasing=F)]
colnames(count)[1:5]

##
umi <- colSums(count)

count.binary <- count 
count.binary[count.binary>0]<-1
geneNum <- colSums(count.binary)

qc$umi <- umi
qc$geneNum <- geneNum


write.table(qc,"RNA.summary.info.new.xls",row.names=F, col.names=T,quote=F,sep="\t")
#summary.info <-data.frame(samplename=name.index,raw_fragment=as.character(df1$V2), clean_fragments=as.character(df2$V2), mapping_rate= df3, uni_mapping_reads=df4, deduplicated_reads=df5)
#summary.info$duplication_rate <-  (summary.info$uni_mapping_reads-summary.info$deduplicated_reads)/summary.info$uni_mapping_reads
