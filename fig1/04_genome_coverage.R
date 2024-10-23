rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F) 

########################################## h3k4me1 (as an example) ###################
setwd("./all_h3k4me1/allBam/05_coverage/")
#cisTopic for countMatrix
library(cisTopic)
library(ggplot2)
library(ggpubr)
library(ggsignif)

pathToBams1 <- c('/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed",invert = T, bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.200bp.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
dat <- apply(c,2,function(x) sum(x>=1))
cover <- (dat/13654391)*100

h3k4me1_cover <- data.frame(cover)
h3k4me1_cover[grep("zygote",row.names(h3k4me1_cover)),2] <- "zygote"
h3k4me1_cover[grep("2cell",row.names(h3k4me1_cover)),2] <- "2cell"
h3k4me1_cover[grep("4cell",row.names(h3k4me1_cover)),2] <- "4cell"
h3k4me1_cover[grep("8cell",row.names(h3k4me1_cover)),2] <- "8cell"
h3k4me1_cover[grep("morula|E2.5",row.names(h3k4me1_cover)),2] <- "morula"
h3k4me1_cover[grep("blastocyst|E3.0|E3.5",row.names(h3k4me1_cover)),2] <- "blastocyst"
colnames(h3k4me1_cover) <- c("coverRate", "stage")

h3k4me1_cover_zygote <- median(h3k4me1_cover[grep("zygote",row.names(h3k4me1_cover)),1])
h3k4me1_cover_2cell <- median(h3k4me1_cover[grep("2cell",row.names(h3k4me1_cover)),1])
h3k4me1_cover_4cell <- median(h3k4me1_cover[grep("4cell",row.names(h3k4me1_cover)),1])
h3k4me1_cover_8cell <- median(h3k4me1_cover[grep("8cell",row.names(h3k4me1_cover)),1])
h3k4me1_cover_morula <- median(h3k4me1_cover[grep("morula",row.names(h3k4me1_cover)),1])
h3k4me1_cover_blastocyst <- median(h3k4me1_cover[grep("blastocyst",row.names(h3k4me1_cover)),1])

h3k4me1_cover_median <- c(h3k4me1_cover_zygote, h3k4me1_cover_2cell, h3k4me1_cover_4cell, h3k4me1_cover_8cell, h3k4me1_cover_morula, h3k4me1_cover_blastocyst)
label <- c("zygote", "2cell", "4cell", "8cell", "morula","blastocyst" )
h3k4me1_cover_median <- data.frame(label, h3k4me1_cover_median)

write.table(h3k4me1_cover, "h3k4me1_coverRate.txt", quote = F)
write.table(h3k4me1_cover_median, "h3k4me1_cover_median.txt", quote = F)
