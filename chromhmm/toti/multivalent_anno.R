########################### multivalent analysis for supple 9d ################################
setwd("./06_multivalent")

library(openxlsx)
########## annotation ############
###annotation with Homer
###annotatePeaks.pl  multivalent.bed  mm10  >  multivalent_anno.txt
anno <- read.xlsx("./multivalent_anno.xlsx")
anno$total <- sum(anno$`Total.size.(bp)`)
anno$ratio <- anno$`Total.size.(bp)`/anno$total

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

pdf("multivalent_anno.pdf")
bp<- ggplot(anno, aes(x="", y=ratio, fill=Annotation))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
bp
dev.off()


######################### alluvial (for supple 9e) #########################

library(networkD3)
library(webshot)
library(tidyr)

###binned bed file should be sorted at first
A <- read.csv("./zygote_inte_all_multivalent_200.txt", header = F, sep = "\t")
B <- read.csv("./2cell_1_inte_all_multivalent_200.txt", header = F, sep = "\t")
C <- read.csv("./2cell_2_inte_all_multivalent_200.txt", header = F, sep = "\t")

A <- unite(A, "min", c("V1", "V2"), sep = ":", remove = T) 
A <- unite(A, "loci", c("min", "V3"), sep = "-", remove = T) 
B <- unite(B, "min", c("V1", "V2"), sep = ":", remove = T) 
B <- unite(B, "loci", c("min", "V3"), sep = "-", remove = T) 
C <- unite(C, "min", c("V1", "V2"), sep = ":", remove = T) 
C <- unite(C, "loci", c("min", "V3"), sep = "-", remove = T) 

identical(A$loci, B$loci)
identical(B$loci, C$loci)
# four cell type ===========================================================================
dat <- cbind(A,B,C)
dat <- dat[,c(1, 2, 4, 6)]
# colnames(dat) <- c("loci","4cell", "8cell", "morula", "TE")
# write.csv(dat, "link_direction_4cell_TE.csv", quote = F)

colnames(dat) <- c("loci","L1", "L2", "L3")
dat$L1 <- paste0("1", dat$L1)
dat$L2 <- paste0("2", dat$L2)
dat$L3 <- paste0("3", dat$L3)
dat$value <- 1
level1 <- aggregate(dat$value, by=list(dat$L1,dat$L2), sum)
level2 <- aggregate(dat$value, by=list(dat$L2,dat$L3), sum)

links <- rbind(level2, level1)
colnames(links) <- c('source',"target","value")
tail(links)

nodes <- data.frame(name=unique(c(links$source,links$target)))
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

Sankey.p <- sankeyNetwork(Links = links, Nodes = nodes,
                          Source = "IDsource", Target = "IDtarget",
                          Value = "value", NodeID = "name", fontFamily = "Arial",
                          fontSize = 12, nodeWidth = 50, nodePadding = 5,
                          height = 200, width = 500,
                          sinksRight=F)
Sankey.p

saveNetwork(Sankey.p, "./multivalent.html", selfcontained = TRUE)
