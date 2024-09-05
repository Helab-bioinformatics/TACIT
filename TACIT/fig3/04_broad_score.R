rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F) 

setwd("./all_h3k36me3/")

library(cisTopic)
library(ggplot2)
library(ggpubr)
library(ggsignif)

pathToBams <- c('./bam2000/')
bamFiles <- paste(pathToBams, list.files(pathToBams), sep='')
bamFiles <- bamFiles[grep("txt|06_reads|bed|bad|agg",invert = T, bamFiles)]

regions1 <- "./broad/broadsubtractsharp.bed"
cisTopicObject1 <- createcisTopicObjectFromBAM(bamFiles, regions1, paired = T, project.name='h3k36me3')
count1 <- as.matrix(cisTopicObject1@count.matrix)
dat1 <- apply(count1,2,function(x) sum(x>=2)) 
h3k36me3_broad <- data.frame(rep("h3k36me3", length(dat1)), dat1)

regions2 <- "./broad/intersect.bed"
cisTopicObject2 <- createcisTopicObjectFromBAM(bamFiles, regions2, paired = T, project.name='h3k36me3')
count2 <- as.matrix(cisTopicObject2@count.matrix)
dat2 <- apply(count2,2,function(x) sum(x>=2)) 
h3k36me3_broad[,3] <- dat2
h3k36me3_broad[,4] <- h3k36me3_broad[,2]+h3k36me3_broad[,3]
h3k36me3_broad[,5] <- h3k36me3_broad[,2]/h3k36me3_broad[,4]

h3k36me3_broad[,1] <- row.names(h3k36me3_broad)
h3k36me3_broad[grep("zygote",h3k36me3_broad[,1]),6] <- "zygote"
h3k36me3_broad[grep("2cell",h3k36me3_broad[,1]),6] <- "2cell"
h3k36me3_broad[grep("4cell",h3k36me3_broad[,1]),6] <- "4cell"
h3k36me3_broad[grep("E2.0",h3k36me3_broad[,1]),6] <- "8cell"
h3k36me3_broad[grep("8cell",h3k36me3_broad[,1]),6] <- "8cell"
h3k36me3_broad[grep("E2.5",h3k36me3_broad[,1]),6] <- "morula"
h3k36me3_broad[grep("morula",h3k36me3_broad[,1]),6] <- "morula"
h3k36me3_broad[grep("E3.0|E3.5|blastocyst",h3k36me3_broad[,1]),6] <- "blastocyst"

h3k36me3_broad <- h3k36me3_broad[,-1]
colnames(h3k36me3_broad) <- c("broad_reads","inte_reads", "all_reads","broad_score","stage")

pdf("./broad/broad.pdf")
p2 <- ggplot(data=h3k36me3_broad,mapping = aes(x = stage, y =broad_score ))+
  geom_violin(aes(fill = stage), trim = F) + 
  ylim(0,1)+
  geom_boxplot(width = 0.05)+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0","#e8490f" ))+
  theme(legend.position = "none")+
  theme_bw()
p2
dev.off()

write.table(h3k36me3_broad,"./h3k36me3_broad.txt", quote = F, col.names = T, row.names = T)

