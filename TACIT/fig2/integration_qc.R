setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")

h3k4me1_2cell <- read.csv("./h3k4me1/h3k4me1_2cell_pseudo.txt", sep = " ")
h3k4me1_4cell <- read.csv("./h3k4me1/h3k4me1_4cell_pseudo.txt", sep = " ")
h3k4me1_8cell <- read.csv("./h3k4me1/h3k4me1_8cell_pseudo.txt", sep = " ")
h3k4me1_morula <- read.csv("./h3k4me1/h3k4me1_morula_pseudo.txt", sep = " ")
h3k4me1_blastocyst <- read.csv("./h3k4me1/h3k4me1_RNA_inte_single.txt", sep = " ")
h3k4me1 <- rbind(h3k4me1_2cell, h3k4me1_4cell, h3k4me1_8cell, h3k4me1_morula, h3k4me1_blastocyst)
h3k4me1$marker <- "H3K4me1"
h3k4me1 <- h3k4me1[which(h3k4me1$prediction.score.max >0.2),]

h3k4me3_2cell <- read.csv("./h3k4me3/h3k4me3_2cell_pseudo_inte.txt", sep = " ")
h3k4me3_4cell <- read.csv("./h3k4me3/h3k4me3_4cell_pseudo_inte.txt", sep = " ")
h3k4me3_8cell <- read.csv("./h3k4me3/h3k4me3_8cell_pseudo_inte.txt", sep = " ")
h3k4me3_morula <- read.csv("./h3k4me3/h3k4me3_morula_pseudo_inte.txt", sep = " ")
h3k4me3_blastocyst <- read.csv("./h3k4me3/h3k4me3_RNA_inte_single.txt", sep = " ")
h3k4me3 <- rbind(h3k4me3_2cell, h3k4me3_4cell, h3k4me3_8cell, h3k4me3_morula, h3k4me3_blastocyst)
h3k4me3$marker <- "h3k4me3"
h3k4me3 <- h3k4me3[which(h3k4me3$prediction.score.max >0.2),]

h3k27ac_2cell <- read.csv("./h3k27ac/h3k27ac_2cell_pseudo.txt", sep = " ")
h3k27ac_4cell <- read.csv("./h3k27ac/h3k27ac_4cell_pseudo.txt", sep = " ")
h3k27ac_8cell <- read.csv("./h3k27ac/h3k27ac_8cell_pseudo.txt", sep = " ")
h3k27ac_morula <- read.csv("./h3k27ac/h3k27ac_morula_pseudo.txt", sep = " ")
h3k27ac_blastocyst <- read.csv("./h3k27ac/h3k27ac_RNA_inte_single.txt", sep = " ")
h3k27ac <- rbind(h3k27ac_2cell, h3k27ac_4cell, h3k27ac_8cell, h3k27ac_morula, h3k27ac_blastocyst)
h3k27ac$marker <- "h3k27ac"
h3k27ac <- h3k27ac[which(h3k27ac$prediction.score.max >0.2),]

h3k36me3_2cell <- read.csv("./h3k36me3/h3k36me3_2cell_pseudo.txt", sep = " ")
h3k36me3_4cell <- read.csv("./h3k36me3/h3k36me3_4cell_pseudo.txt", sep = " ")
h3k36me3_8cell <- read.csv("./h3k36me3/h3k36me3_8cell_pseudo.txt", sep = " ")
h3k36me3_morula <- read.csv("./h3k36me3/h3k36me3_morula_pseudo.txt", sep = " ")
h3k36me3_blastocyst <- read.csv("./h3k36me3/h3k36me3_RNA_inte_single.txt", sep = " ")
h3k36me3 <- rbind(h3k36me3_2cell, h3k36me3_4cell, h3k36me3_8cell, h3k36me3_morula, h3k36me3_blastocyst)
h3k36me3$marker <- "h3k36me3"
h3k36me3 <- h3k36me3[which(h3k36me3$prediction.score.max >0.2),]

cotacit_2cell <- read.csv("./cotacit_data/cotacit_2cell_h3k27ac_pseudo.txt", sep = " ")
cotacit_4cell <- read.csv("./cotacit_data/cotacit_4cell_h3k27ac_pseudo.txt", sep = " ")
cotacit_8cell <- read.csv("./cotacit_data/cotacit_8cell_h3k27ac_pseudo.txt", sep = " ")
cotacit_morula <- read.csv("./cotacit_data/cotacit_morula_h3k27ac_pseudo.txt", sep = " ")
cotacit_blastocyst <- read.csv("./cotacit_data/cotacit_blasto_h3k27ac_RNA_inte_single.txt", sep = " ")
cotacit <- rbind(cotacit_2cell, cotacit_4cell, cotacit_8cell, cotacit_morula, cotacit_blastocyst)
cotacit$marker <- "coTACIT"
cotacit <- cotacit[which(cotacit$prediction.score.max >0.2),]

h3k36me3_8cell$marker <- "H3K36me3"
h3k4me3_8cell$marker <- "H3K4me3"
h3k4me1_8cell$marker <- "H3K4me1"
h3k27ac_8cell$marker <- "H3K27ac"
cotacit_8cell$marker <- "coTACIT"

#all_inte <- rbind(h3k4me1, h3k4me3, h3k27ac, h3k36me3, cotacit)
all_inte <- rbind(h3k4me1_8cell, h3k4me3_8cell, h3k27ac_8cell, h3k36me3_8cell, cotacit_8cell)
colnames(all_inte) <- c("predicted_id", "score", "marker")
all_inte <- all_inte[which(all_inte$score >0.2),]
all_inte$marker <- factor(all_inte$marker, levels = c("H3K4me1", "H3K4me3", "H3K27ac", "H3K36me3", "coTACIT"))

library(ggplot2)
pdf("./all_marker_8cell_integration_score.pdf")
ggplot(data=all_inte,mapping = aes(x = marker, y = score))+
  geom_violin(aes(fill = marker), trim = FALSE) + 
  geom_jitter()+
  ylim(0,1)+
  #geom_boxplot(width = 0.05)+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0"))+
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
dev.off()

