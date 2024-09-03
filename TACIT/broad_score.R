rm(list = ls())  ## 魔幻操作，一键清空~
Sys.setenv(R_MAX_NUM_DLLS=999) ##Sys.setenv修改环境设置，R的namespace是有上限的，如果导入包时超过这个上次就会报错,R_MAX_NUM_DLLS可以修改这个上限
options(stringsAsFactors = F) ##options:允许用户对工作空间进行全局设置，stringsAsFactors防止R自动把字符串string的列辨认成factor

# options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

setwd("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/")

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


p2 <- ggplot(data=h3k36me3_broad,mapping = aes(x = stage, y =broad_score ))+
  geom_violin(aes(fill = stage), trim = F) + 
  ylim(0,1)+
  geom_boxplot(width = 0.05)+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0","#e8490f" ))+
  theme(legend.position = "none")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
pdf("./broad/broad.pdf")
p2
dev.off()

write.table(h3k36me3_broad,"./h3k36me3_broad.txt", quote = F, col.names = T, row.names = T)

