#setwd("D:/Data_Min/tacit/figures/24_allstage_integration/blasto")
setwd("D:/Data_Min/tacit/figures/24_allstage_integration/blasto/randomForest/top780")
library(tidyr)
library(ggplot2)
library(openxlsx)


######################################### ICM ##########################################
#setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/learnmodel_2kb/POSTERIOR/test/states/ICM")

library(tidyr)
library(ggplot2)
ICM <- read.csv("./ICM_enhancer_states_motifs.csv", sep = ",")
ICM_4cell <- ICM[grep("4cell", ICM$stage),]
ICM_8cell <- ICM[grep("8cell", ICM$stage),]
ICM_morula <- ICM[grep("morula", ICM$stage),]
ICM_blasto <- ICM[grep("ICM", ICM$stage),]

ICM_4cell <- ICM_4cell[!duplicated(ICM_4cell$TF),]
ICM_8cell <- ICM_8cell[!duplicated(ICM_8cell$TF),]
ICM_morula <- ICM_morula[!duplicated(ICM_morula$TF),]
ICM_blasto <- ICM_blasto[!duplicated(ICM_blasto$TF),]
ICM <- rbind(ICM_4cell, ICM_8cell, ICM_morula, ICM_blasto)

# ICM_inte_all <- Reduce(intersect, list(ICM_4cell$TF, ICM_8cell$TF, ICM_morula$TF, ICM_blasto$TF))
# ICM_inte_48 <- Reduce(intersect, list(ICM_4cell$TF, ICM_8cell$TF))
# ICM_inte_8m <- Reduce(intersect, list(ICM_8cell$TF, ICM_morula$TF))
# ICM_inte_mb <- Reduce(intersect, list(ICM_morula$TF, ICM_blasto$TF))
# ICM_inte_no4 <- Reduce(intersect, list(ICM_8cell$TF, ICM_morula$TF, ICM_blasto$TF))
# 
# ICM_inte_all <- ICM[which(ICM$TF %in% ICM_inte_all),]
# ICM_inte_48 <- ICM[which(ICM$TF %in% ICM_inte_48),]
# ICM_inte_8m <- ICM[which(ICM$TF %in% ICM_inte_8m),]
# ICM_inte_mb <- ICM[which(ICM$TF %in% ICM_inte_mb),]
# ICM_inte_no4 <- ICM[which(ICM$TF %in% ICM_inte_no4),]


select_motifs <- c("DUX",  "NR5A2", "ZSCAN4D", 
                   "POU5F1", "NANOG",  "FOXA2","GATA4",
                   "ZFX","HNF4A","YY2","TCF12","CEBPB", "BBX", "SMAD2","HBP1","PRDM14", "PRDM15","CDX2","KLF6","SOX15","MED1","ELF5","HIF1A","KLF5","ELF3", "EWSR1")
selected <- ICM[which(ICM$TF %in% select_motifs),]

genes <- read.csv("D:/Data_Min/tacit/reference/genes/2014_science_RNA/GSE45719_RPKM_scRNA.csv", sep = ",")
genes <- genes[!duplicated(genes$Gene_symbol),]
rownames(genes) <- genes[,1]
genes <- genes[,-1]
rownames(genes) <- toupper(rownames(genes))

ICM_id <- read.csv("./blasto_pseudo_ICM.txt", sep = " ")
TE_id <- read.csv("./blasto_pseudo_TE.txt", sep = " ")
a <- intersect(rownames(ICM_id), rownames(TE_id))

genes_4cell <- genes[,grep("X4cell_2.3|X4cell_3.1|X4cell_3.3|X4cell_2.1|X4cell_2.2|X4cell_2.4|X4cell_3.4", colnames(genes))]
genes_8cell <- genes[,grep("X8cell_1.7|X8cell_1.8|X8cell_1.5|X8cell_1.6|X8cell_1.1|X8cell_1.2|X8cell_1.4|X8cell_8.3|X8cell_8.2|X8cell_8.7", colnames(genes))]
genes_morula <- genes[,grep("X16cell_5.6|X16cell_6.3|X16cell_6.4|X16cell_6.6|X16cell_6.7|X16cell_6.8|X16cell_6.12|X16cell_6.1|X16cell_6.2|X16cell_6.11", invert = T,colnames(genes))]
genes_ICM <- genes[,which(colnames(genes) %in% rownames(ICM_id))]

exp <- data.frame(genes_4cell=rowMeans(genes_4cell),
                  genes_8cell=rowMeans(genes_8cell),
                  genes_morula=rowMeans(genes_morula),
                  genes_ICM=rowMeans(genes_ICM))
exp <- exp[which(rownames(exp) %in% select_motifs),]
exp$gene <- rownames(exp)
exp$TF <- rownames(exp)
exp <- gather(exp, key=cell, value=exp, -c(gene, TF))
exp <- unite(exp, id, c("cell", "TF"), sep = "_", remove = F)

selected <- selected[grep("DCE", invert = T, selected$TF),]
selected$genes <- "genes"
selected <- unite(selected, id, c("genes", "stage", "TF"), sep = "_", remove = F)
exp$pvalue <- selected$pvalue[match(exp$id, selected$id)]

exp$TF <- factor(exp$TF, levels = select_motifs)
exp$cell <- factor(exp$cell, levels = c("genes_4cell", "genes_8cell", "genes_morula", "genes_ICM"))

exp[is.na(exp)] <- 0
exp[which(exp$exp == 0),7] <- "L1"
exp[which(exp$exp >0 & exp$exp<10),7] <- "L2"
exp[which(exp$exp >10 & exp$exp<20),7] <- "L3"
exp[which(exp$exp >20 & exp$exp<30),7] <- "L4"
exp[which(exp$exp >30 ),7] <- "L5"
colnames(exp) <- c(colnames(exp[,1:6]), "group")

exp$logexp <- log2(exp$exp + 1)
#exp[which(exp$logexp > 1.5),8] <- 1.5

exp[which(exp$logexp == 0),9] <- "L1"
exp[which(exp$logexp >0 & exp$logexp<1),9] <- "L2"
exp[which(exp$logexp >1 & exp$logexp<3),9] <- "L3"
exp[which(exp$logexp >3 & exp$logexp<5.1),9] <- "L4"
exp[which(exp$logexp >5.1 ),9] <- "L5"
colnames(exp) <- c(colnames(exp[,1:8]), "level")

exp[is.na(exp)] <- 0
pdf("ICM_motifs_dot.pdf", width = 3, height = 5)
ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)

ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`level`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()

write.table(exp, "ICM_selected_TFs.txt", quote = F, sep = "\t")

######################################### TE ##########################################
#setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/learnmodel_2kb/POSTERIOR/test/states/TE")

library(tidyr)
library(ggplot2)
TE <- read.csv("./TE_enhancer_states_motifs.csv", sep = ",")
TE_4cell <- TE[grep("4cell", TE$stage),]
TE_8cell <- TE[grep("8cell", TE$stage),]
TE_morula <- TE[grep("morula", TE$stage),]
TE_blasto <- TE[grep("TE", TE$stage),]

TE_4cell <- TE_4cell[!duplicated(TE_4cell$TF),]
TE_8cell <- TE_8cell[!duplicated(TE_8cell$TF),]
TE_morula <- TE_morula[!duplicated(TE_morula$TF),]
TE_blasto <- TE_blasto[!duplicated(TE_blasto$TF),]
TE <- rbind(TE_4cell, TE_8cell, TE_morula, TE_blasto)

# genes <- read.csv("D:/Data_Min/tacit/reference/genes/2014_science_RNA/GSE45719_RPKM_scRNA.csv", sep = ",")
# genes <- genes[!duplicated(genes$Gene_symbol),]
# rownames(genes) <- genes[,1]
# genes <- genes[,-1]
# genes_TE <- genes[,which(colnames(genes) %in% rownames(TE_id))]
# genes_TE <- apply(genes_TE, 1, median)
# TE_blasto$exp <- genes_TE[match(tolower(TE_blasto$TF), tolower(names(genes_TE)))]
# TE_select <- TE_blasto[which(TE_blasto$pvalue > 5 & TE_blasto$exp >8),]
# TE_2 <- TE[which(TE$TF %in% TE_select$TF),]
# 
# genes_ICM <- genes[,which(colnames(genes) %in% rownames(ICM_id))]
# genes_ICM <- apply(genes_ICM, 1, median)
# ICM_blasto$exp <- genes_ICM[match(tolower(ICM_blasto$TF), tolower(names(genes_ICM)))]
# ICM_select <- ICM_blasto[which(ICM_blasto$pvalue > 5 & ICM_blasto$exp >8),]



# TE_inte_all <- Reduce(intersect, list(TE_4cell$TF, TE_8cell$TF, TE_morula$TF, TE_blasto$TF))
# TE_inte_48 <- Reduce(intersect, list(TE_4cell$TF, TE_8cell$TF))
# TE_inte_8m <- Reduce(intersect, list(TE_8cell$TF, TE_morula$TF))
# TE_inte_mb <- Reduce(intersect, list(TE_morula$TF, TE_blasto$TF))
# TE_inte_no4 <- Reduce(intersect, list(TE_8cell$TF, TE_morula$TF, TE_blasto$TF))
# 
# TE_inte_all <- TE[which(TE$TF %in% TE_inte_all),]
# TE_inte_48 <- TE[which(TE$TF %in% TE_inte_48),]
# TE_inte_8m <- TE[which(TE$TF %in% TE_inte_8m),]
# TE_inte_mb <- TE[which(TE$TF %in% TE_inte_mb),]
# TE_inte_no4 <- TE[which(TE$TF %in% TE_inte_no4),]

selected <- TE[which(TE$TF %in% select_motifs),]

genes_4cell <- genes[,grep("X4cell_4.2|X4cell_4.3|X4cell_1.4", colnames(genes))]
genes_8cell <- genes[,grep("X8cell_5.4|X8cell_5.1|X8cell_5.2|X8cell_5.3|X8cell_5.8|X8cell_5.6|X8cell_5.7|X8cell_8.8|X8cell_8.1|X8cell_8.4|X8cell_8.6", colnames(genes))]
genes_morula <- genes[,grep("X16cell_5.6|X16cell_6.3|X16cell_6.4|X16cell_6.6|X16cell_6.7|X16cell_6.8|X16cell_6.12|X16cell_6.1|X16cell_6.2|X16cell_6.11", colnames(genes))]
genes_TE <- genes[,which(colnames(genes) %in% rownames(TE_id))]

exp <- data.frame(genes_4cell=apply(genes_4cell, 1, median),
                   genes_8cell=apply(genes_8cell, 1, median),
                   genes_morula=apply(genes_morula,1, median),
                   genes_TE=apply(genes_TE, 1, median))
exp <- exp[which(rownames(exp) %in% select_motifs),]
exp$gene <- rownames(exp)
exp$TF <- exp$gene
exp <- gather(exp, key=cell, value=exp, -c(gene, TF))
exp <- unite(exp, id, c("cell", "TF"), sep = "_", remove = F)

selected <- selected[grep("DCE", invert = T, selected$TF),]
selected$genes <- "genes"
selected <- unite(selected, id, c("genes", "stage", "TF"), sep = "_", remove = F)
exp$pvalue <- selected$pvalue[match(exp$id, selected$id)]

exp$TF <- factor(exp$TF, levels = select_motifs)
exp$cell <- factor(exp$cell, levels = c("genes_4cell", "genes_8cell", "genes_morula", "genes_TE"))

exp[which(exp$exp == 0),7] <- "L1"
exp[which(exp$exp >0 & exp$exp<10),7] <- "L2"
exp[which(exp$exp >10 & exp$exp<20),7] <- "L3"
exp[which(exp$exp >20 & exp$exp<30),7] <- "L4"
exp[which(exp$exp >30 ),7] <- "L5"
colnames(exp) <- c(colnames(exp[,1:6]), "group")

exp$logexp <- log2(exp$exp + 1)
#exp[which(exp$logexp > 1.5),8] <- 1.5

exp[which(exp$logexp == 0),9] <- "L1"
exp[which(exp$logexp >0 & exp$logexp<1),9] <- "L2"
exp[which(exp$logexp >1 & exp$logexp<3),9] <- "L3"
exp[which(exp$logexp >3 & exp$logexp<5.1),9] <- "L4"
exp[which(exp$logexp >5.1 ),9] <- "L5"
colnames(exp) <- c(colnames(exp[,1:8]), "level")

exp[is.na(exp)] <- 0
pdf("TE_motifs_dot.pdf", width = 3, height = 5)
ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)

ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`level`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()

write.table(exp, "TE_selected_TFs.txt", quote = F, sep = "\t")




############################## select more motifs for screening ###########################
library(tidyr)
library(ggplot2)
ICM <- read.csv("./ICM_enhancer_states_motifs.csv", sep = ",")
ICM_4cell <- ICM[grep("4cell", ICM$stage),]
ICM_8cell <- ICM[grep("8cell", ICM$stage),]
ICM_morula <- ICM[grep("morula", ICM$stage),]
ICM_blasto <- ICM[grep("ICM", ICM$stage),]

ICM_4cell <- ICM_4cell[!duplicated(ICM_4cell$TF),]
ICM_8cell <- ICM_8cell[!duplicated(ICM_8cell$TF),]
ICM_morula <- ICM_morula[!duplicated(ICM_morula$TF),]
ICM_blasto <- ICM_blasto[!duplicated(ICM_blasto$TF),]

TE <- read.csv("./TE_enhancer_states_motifs.csv", sep = ",")
TE_4cell <- TE[grep("4cell", TE$stage),]
TE_8cell <- TE[grep("8cell", TE$stage),]
TE_morula <- TE[grep("morula", TE$stage),]
TE_blasto <- TE[grep("TE", TE$stage),]

TE_4cell <- TE_4cell[!duplicated(TE_4cell$TF),]
TE_8cell <- TE_8cell[!duplicated(TE_8cell$TF),]
TE_morula <- TE_morula[!duplicated(TE_morula$TF),]
TE_blasto <- TE_blasto[!duplicated(TE_blasto$TF),]

##### for ICM #####
TFs_ICM_morula <- ICM_morula[which(ICM_morula$pvalue >2),"TF"]
TFs_ICM_blasto <- ICM_blasto[which(ICM_blasto$pvalue >2),"TF"]
TFs_ICM <- intersect(TFs_ICM_morula, TFs_ICM_blasto)

no_blasto <- TE_blasto$pvalue[match(TFs_ICM, TE_blasto$TF)]
names(no_blasto) <- TFs_ICM
no_blasto[is.na(no_blasto)] <- 0
TFs_ICM <- names(no_blasto)[no_blasto < 3]

##### for TE #####
TFs_TE_morula <- TE_morula[which(TE_morula$pvalue >2),"TF"]
TFs_TE_blasto <- TE_blasto[which(TE_blasto$pvalue >2),"TF"]
TFs_TE <- intersect(TFs_TE_morula, TFs_TE_blasto)

no_blasto <- ICM_blasto$pvalue[match(TFs_TE, ICM_blasto$TF)]
names(no_blasto) <- TFs_TE
no_blasto[is.na(no_blasto)] <- 0
TFs_TE <- names(no_blasto)[no_blasto < 3]

##### combine #####
TFs <- union(TFs_ICM, TFs_TE)

TFs_ICM_4cell <- ICM_4cell$pvalue[match(TFs, ICM_4cell$TF)]
TFs_ICM_8cell <- ICM_8cell$pvalue[match(TFs, ICM_8cell$TF)]
TFs_ICM_morula <- ICM_morula$pvalue[match(TFs, ICM_morula$TF)]
TFs_ICM_blasto <- ICM_blasto$pvalue[match(TFs, ICM_blasto$TF)]

TFs_TE_4cell <- TE_4cell$pvalue[match(TFs, TE_4cell$TF)]
TFs_TE_8cell <- TE_8cell$pvalue[match(TFs, TE_8cell$TF)]
TFs_TE_morula <- TE_morula$pvalue[match(TFs, TE_morula$TF)]
TFs_TE_blasto <- TE_blasto$pvalue[match(TFs, TE_blasto$TF)]

TFs_ICM_4cell[is.na(TFs_ICM_4cell)] <- 0
TFs_ICM_8cell[is.na(TFs_ICM_8cell)] <- 0
TFs_ICM_morula[is.na(TFs_ICM_morula)] <- 0
TFs_ICM_blasto[is.na(TFs_ICM_blasto)] <- 0
TFs_TE_4cell[is.na(TFs_TE_4cell)] <- 0
TFs_TE_8cell[is.na(TFs_TE_8cell)] <- 0
TFs_TE_morula[is.na(TFs_TE_morula)] <- 0
TFs_TE_blasto[is.na(TFs_TE_blasto)] <- 0

dat <- data.frame(ICM_4cell=TFs_ICM_4cell, ICM_8cell=TFs_ICM_8cell, ICM_morula=TFs_ICM_morula, ICM_blasto=TFs_ICM_blasto, 
                  TE_4cell=TFs_TE_4cell, TE_8cell=TFs_TE_8cell, TE_morula=TFs_TE_morula, TE_blasto=TFs_TE_blasto)
rownames(dat) <- TFs

write.table(dat, "ICM_TE_diff_TFs_list.txt", sep =" ", quote = F)

############################### pairing ms, Ribo-seq and RNA-seq data ##################################
setwd("D:/Data_Min/tacit/figures/24_allstage_integration/blasto/randomForest/top780")
#install.packages("readxl")
library(readxl)

tfs <- read.csv("ICM_TE_diff_TFs_list.txt", sep = " ")


ribo <- read_excel("D:/Data_Min/tacit/reference/TFs/2022_Ribo-lite_xiewei.xlsx", sheet = 1)
ribo <- ribo[,1:12]
colnames(ribo) <- ribo[1,]
ribo <- ribo[-1,]
tfs$C4_ribo <- ribo$`4C`[match(tolower(rownames(tfs)), tolower(ribo$gene))]
tfs$C8_ribo <- ribo$`8C`[match(tolower(rownames(tfs)), tolower(ribo$gene))]
tfs$ICM_ribo <- ribo$ICM[match(tolower(rownames(tfs)), tolower(ribo$gene))]

ms <- read_excel("D:/Data_Min/tacit/reference/TFs/2023_ms.xlsx")
tfs$ms_C4_rep1 <- ms$`4C_rep1`[match(tolower(rownames(tfs)), tolower(ms$`gene name`))]
tfs$ms_C4_rep2 <- ms$`4C_rep2`[match(tolower(rownames(tfs)), tolower(ms$`gene name`))]
tfs$ms_C8_rep1 <- ms$`8C_rep1`[match(tolower(rownames(tfs)), tolower(ms$`gene name`))]
tfs$ms_C8_rep2 <- ms$`8C_rep2`[match(tolower(rownames(tfs)), tolower(ms$`gene name`))]
tfs$ms_blasto_rep1 <- ms$Blastocyst_rep1[match(tolower(rownames(tfs)), tolower(ms$`gene name`))]
tfs$ms_blasto_rep2 <- ms$Blastocyst_rep2[match(tolower(rownames(tfs)), tolower(ms$`gene name`))]


genes <- read.csv("D:/Data_Min/tacit/reference/genes/2014_science_RNA/GSE45719_RPKM_scRNA.csv", sep = ",")
genes <- genes[!duplicated(genes$Gene_symbol),]
rownames(genes) <- genes[,1]
genes <- genes[,-1]
selected_rows <- tolower(rownames(genes)) %in% tolower(rownames(tfs))
genes_selected <- genes[selected_rows, ]

genes_4cell <- genes_selected[,grep("4cell", colnames(genes_selected))]
genes_8cell <- genes_selected[,grep("8cell", colnames(genes_selected))]
genes_morula <- genes_selected[,grep("16cell", colnames(genes_selected))]
genes_blasto <- genes_selected[,grep("blast", colnames(genes_selected))]

genes_4cell <- rowMeans(genes_4cell)
genes_8cell <- rowMeans(genes_8cell)
genes_morula <- rowMeans(genes_morula)
genes_blasto <- rowMeans(genes_blasto)

exp <- data.frame(genes_4cell, genes_8cell, genes_morula, genes_blasto)

tfs[,18:21] <- exp[match(tolower(rownames(tfs)), tolower(rownames(exp))),1:4]

tfs <- tfs[,c(1:8, 18:21, 9:17)]
write.table(tfs, "ICM_TE_diff_TFs_list_ribo.txt", sep =" ", quote = F)


######################################## 2024.2.7 ####################################



