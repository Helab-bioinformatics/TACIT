
######################################### ICM ##########################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/learnmodel_2kb/POSTERIOR/test/states/ICM")

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

inte_all <- Reduce(intersect, list(ICM_4cell$TF, ICM_8cell$TF, ICM_morula$TF, ICM_blasto$TF))
inte_48 <- Reduce(intersect, list(ICM_4cell$TF, ICM_8cell$TF))
inte_8m <- Reduce(intersect, list(ICM_8cell$TF, ICM_morula$TF))
inte_mb <- Reduce(intersect, list(ICM_morula$TF, ICM_blasto$TF))
inte_no4 <- Reduce(intersect, list(ICM_8cell$TF, ICM_morula$TF, ICM_blasto$TF))

inte_all <- ICM[which(ICM$TF %in% inte_all),]
inte_48 <- ICM[which(ICM$TF %in% inte_48),]
inte_8m <- ICM[which(ICM$TF %in% inte_8m),]
inte_mb <- ICM[which(ICM$TF %in% inte_mb),]
inte_no4 <- ICM[which(ICM$TF %in% inte_no4),]


select_motifs <- c("DUXBL", "NR5A1", "NR5A2", "ZSCAN4", "RARG", 
                    "SOX2", "SOX3","POU5F1", "NANOG", "YY2", "GATA1", "GATA2", "GATA6", "EBF2", "EBF3",
                   "SMAD3","DCE", "TCF3","CEBPB")
selected <- ICM[which(ICM$TF %in% select_motifs),]

genes <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/GSE45719_scRNA.csv", sep = ",")
genes <- genes[grep("Duxbl|Nr5a1|Nr5a2|Zscan4|Rarg|Sox2|Sox3|Pou5f1|Nanog|Yy2|Gata1|Gata2|Gata6|Smad3|Tcf3|Cebpb|Ebf2|Ebf3", genes$Gene_symbol),]
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- genes[grep("Nanogpd|Sox21|Sox2ot|Sox30|Zscan4a|Zscan4b|Zscan4d|Zscan4f|Zscan4e|Gm10391", invert = T, rownames(genes)),]

ICM_id <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/blasto_pseudo_ICM.txt", sep = " ")
TE_id <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/blasto_pseudo_TE.txt", sep = " ")
a <- intersect(rownames(ICM_id), rownames(TE_id))

genes_4cell <- genes[,grep("4cell", colnames(genes))]
genes_8cell <- genes[,grep("8cell", colnames(genes))]
genes_morula <- genes[,grep("16cell", colnames(genes))]
genes_ICM <- genes[,which(colnames(genes) %in% rownames(ICM_id))]

exp <- data.frame(genes_4cell=rowMeans(genes_4cell),
                  genes_8cell=rowMeans(genes_8cell),
                  genes_morula=rowMeans(genes_morula),
                  genes_ICM=rowMeans(genes_ICM))
exp$gene <- rownames(exp)
exp$TF <- c("TCF3", "NR5A2", "EBF3","NANOG", "ZSCAN4", "SOX3", "EBF2","POU5F1", "RARG",
            "DUXBL", "GATA6", "SMAD3", "CEBPB", "GATA2", "GATA1", "NR5A1", "SOX2", "YY2")
exp <- gather(exp, key=cell, value=exp, -c(gene, TF))
exp <- unite(exp, id, c("cell", "TF"), sep = "_", remove = F)

selected <- selected[grep("DCE", invert = T, selected$TF),]
selected$genes <- "genes"
selected <- unite(selected, id, c("genes", "stage", "TF"), sep = "_", remove = F)
exp$pvalue <- selected$pvalue[match(exp$id, selected$id)]

exp$TF <- factor(exp$TF, levels = c("SMAD3", "TCF3","CEBPB",  
                                    "SOX3", 'EBF2',"EBF3","YY2", "GATA1", "GATA2", "GATA6", "SOX2", "POU5F1", "NANOG",
                                    "DUXBL", "NR5A1", "NR5A2", "ZSCAN4", "RARG"))
exp$cell <- factor(exp$cell, levels = c("genes_4cell", "genes_8cell", "genes_morula", "genes_ICM"))

exp[which(exp$exp == 0),7] <- "L1"
exp[which(exp$exp >0 & exp$exp<5),7] <- "L2"
exp[which(exp$exp >5 & exp$exp<10),7] <- "L3"
exp[which(exp$exp >10 & exp$exp<20),7] <- "L4"
exp[which(exp$exp >20 ),7] <- "L5"
colnames(exp) <- c(colnames(exp[,1:6]), "group")

exp$logexp <- log10(exp$exp + 1)
exp[which(exp$logexp > 1.5),8] <- 1.5

pdf("ICM_motifs_dot.pdf", width = 4, height = 5)
ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)

ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`logexp`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()

pdf("ICM_motifs_dot_2.pdf", width = 10, height = 3)
ggplot(exp,aes(x=TF,y=cell))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()


######################################### TE ##########################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/learnmodel_2kb/POSTERIOR/test/states/TE")

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

inte_all <- Reduce(intersect, list(TE_4cell$TF, TE_8cell$TF, TE_morula$TF, TE_blasto$TF))
inte_48 <- Reduce(intersect, list(TE_4cell$TF, TE_8cell$TF))
inte_8m <- Reduce(intersect, list(TE_8cell$TF, TE_morula$TF))
inte_mb <- Reduce(intersect, list(TE_morula$TF, TE_blasto$TF))
inte_no4 <- Reduce(intersect, list(TE_8cell$TF, TE_morula$TF, TE_blasto$TF))

inte_all <- TE[which(TE$TF %in% inte_all),]
inte_48 <- TE[which(TE$TF %in% inte_48),]
inte_8m <- TE[which(TE$TF %in% inte_8m),]
inte_mb <- TE[which(TE$TF %in% inte_mb),]
inte_no4 <- TE[which(TE$TF %in% inte_no4),]

select_motifs <- c("TCF3", "SOX17","OSR1","ZIC1",
                   "FOXK1","HMBOX1","IRF6","SOX13","SOX15", "KLF6","MEIS2", "NFIA", "TFAP2A", "TFAP2C","CDX2","TEAD4",
                   "DUX","NR5A1", "NR5A2")
selected <- TE[which(TE$TF %in% select_motifs),]

genes <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/GSE45719_scRNA.csv", sep = ",")
genes <- genes[grep("Dux|Nr5a1|Nr5a2|Cdx2|Tead4|Klf6|Meis2|Sox17|Nfia|Tfap2a|Tfap2c|Tcf3|Osr1|Zic1|Foxk1|Hmbox1|Irf6|Sox13|Sox15", genes$Gene_symbol),]
rownames(genes) <- genes[,1]
genes <- genes[,-1]
genes <- genes[grep("Duxbl", invert = T, rownames(genes)),]


TE_id <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/blasto_pseudo_TE.txt", sep = " ")

genes_4cell <- genes[,grep("4cell", colnames(genes))]
genes_8cell <- genes[,grep("8cell", colnames(genes))]
genes_morula <- genes[,grep("16cell", colnames(genes))]
genes_TE <- genes[,which(colnames(genes) %in% rownames(TE_id))]

exp <- data.frame(genes_4cell=rowMeans(genes_4cell),
                  genes_8cell=rowMeans(genes_8cell),
                  genes_morula=rowMeans(genes_morula),
                  genes_TE=rowMeans(genes_TE))
exp$gene <- rownames(exp)
exp$TF <- c("TCF3", "ZIC1","TEAD4", "NR5A2","SOX15", "HMBOX1","IRF6","KLF6", "NFIA", "MEIS2", "DUX", "SOX13","CDX2", "NR5A1", "FOXK1","OSR1", "SOX17")
exp <- gather(exp, key=cell, value=exp, -c(gene, TF))
exp <- unite(exp, id, c("cell", "TF"), sep = "_", remove = F)

selected <- selected[grep("DCE", invert = T, selected$TF),]
selected$genes <- "genes"
selected <- unite(selected, id, c("genes", "stage", "TF"), sep = "_", remove = F)
exp$pvalue <- selected$pvalue[match(exp$id, selected$id)]

exp$TF <- factor(exp$TF, levels = c("TCF3", "SOX17","OSR1","ZIC1",
                                     "FOXK1","HMBOX1","IRF6","SOX13","SOX15", "KLF6","MEIS2", "NFIA", "TFAP2A", "TFAP2C","CDX2","TEAD4",
                                    "DUX","NR5A1", "NR5A2"))
exp$cell <- factor(exp$cell, levels = c("genes_4cell", "genes_8cell", "genes_morula", "genes_TE"))

exp[which(exp$exp == 0),7] <- "L1"
exp[which(exp$exp >0 & exp$exp<5),7] <- "L2"
exp[which(exp$exp >5 & exp$exp<10),7] <- "L3"
exp[which(exp$exp >10 & exp$exp<20),7] <- "L4"
exp[which(exp$exp >20 ),7] <- "L5"
colnames(exp) <- c(colnames(exp[,1:6]), "group")

exp$logexp <- log10(exp$exp + 1)
exp[which(exp$logexp > 1.5),8] <- 1.5

pdf("TE_motifs_dot.pdf", width = 4, height = 5)
ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)

ggplot(exp,aes(x=cell,y=TF))+
  geom_point(aes(size=`pvalue`,
                 color=`logexp`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()

pdf("TE_motifs_dot_2.pdf", width = 10, height = 3)
ggplot(exp,aes(x=TF,y=cell))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()

