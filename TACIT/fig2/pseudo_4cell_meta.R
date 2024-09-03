##pseudo_4cell_integration.R脚本把map到同一个RNA pseudocells的多个TACIT细胞只取了一个，这导致信号稀疏，track图不好看；
##本脚本把map到同一个RNA pseudocells的多个TACIT细胞进行aggregate，合成了meta cells；所以整合的时候实际上是meta cells；这样保留的信号更多；

library(tidyr)
####################################### all #######################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_4cell_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_4cell_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_4cell_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_4cell_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_4cell_h3k27ac_pseudo.txt", sep = " ")

h3k4me1_pre1 <- h3k4me1_pre1[which(h3k4me1_pre1$prediction.score.max >0.2),]
h3k4me3_pre1 <- h3k4me3_pre1[which(h3k4me3_pre1$prediction.score.max >0.2),]
h3k27ac_pre1 <- h3k27ac_pre1[which(h3k27ac_pre1$prediction.score.max >0.2),]
h3k36me3_pre1 <- h3k36me3_pre1[which(h3k36me3_pre1$prediction.score.max >0.2),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))

#### h3k4me1 ####
dir.sh <- paste0("mkdir ", inte)
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1 $predicted.id %in% inte),]
sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/", rownames(h3k4me1_pre1), 
                     " ./", h3k4me1_pre1$predicted.id, "/")
h3k4me1.sh <- c(dir.sh, sh)
write.table(h3k4me1.sh, "./03_chromHMM_4cell_meta/h3k4me1_4cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k4me3 ####
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", rownames(h3k4me3_pre1), 
               " ./", h3k4me3_pre1$predicted.id, "/")
write.table(h3k4me3.sh, "./03_chromHMM_4cell_meta/h3k4me3_4cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k27ac ####
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", rownames(h3k27ac_pre1), 
               " ./", h3k27ac_pre1$predicted.id, "/")
write.table(h3k27ac.sh, "./03_chromHMM_4cell_meta/h3k27ac_4cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k36me3 ####
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", rownames(h3k36me3_pre1), 
               " ./", h3k36me3_pre1$predicted.id, "/")
write.table(h3k36me3.sh, "./03_chromHMM_4cell_meta/h3k36me3_4cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)


#### cotacit ####
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_4cell_h3k27ac_pseudo.txt", sep = " ")
cotacit_h3k27me3_pre1 <- read.csv("./cotacit_data/h3k27me3_4cell_cotacit_tacit_inte.txt", sep = " ")
cotacit_h3k9me3_pre1 <- read.csv("./cotacit_data/h3k9me3_4cell_cotacit_tacit_inte.txt", sep = " ")

cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1 $predicted.id %in% inte),]
cotacit_pre1$ID <- rownames(cotacit_pre1)
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_all_h3k9me3_rmdup.bam") 
cotacit_pre1$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1$ID3, rownames(cotacit_h3k9me3_pre1))]

k27me3_tacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1$k27_tacit, 
                      " ./", cotacit_pre1$predicted.id, "/")
k27me3_cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2, 
                          " ./", cotacit_pre1$predicted.id, "/")
k9me3_tacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1$k9_tacit, 
                          " ./", cotacit_pre1$predicted.id, "/")
k9me3_cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3, 
                            " ./", cotacit_pre1$predicted.id, "/")
k27me3.sh <- c(k27me3_tacit.sh, k27me3_cotacit.sh)
k9me3.sh <- c(k9me3_tacit.sh, k9me3_cotacit.sh)
write.table(k27me3.sh, "./03_chromHMM_4cell_meta/h3k27me3_4cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./03_chromHMM_4cell_meta/h3k9me3_4cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

write.table(cotacit_pre1, "./03_chromHMM_4cell_meta/cotacit_tacit_inte_4cell.txt", quote = F, sep = " ")
