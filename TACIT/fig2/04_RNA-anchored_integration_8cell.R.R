
#######################################################################################################################
####################################### pseudo-cell integration #######################################################
library(cisTopic)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(edgeR)
library(openxlsx)
set.seed(1)

####################################### 8cell integration ######################################################
RNA_pseudo <- readRDS("./scRNA/RNA_pseudo.rds")
RNA_pseudo <- subset(RNA_pseudo,cells = colnames(RNA_pseudo[,grep("8cell", RNA_pseudo@meta.data$group.ident)]))
RNA_pseudo <- NormalizeData(RNA_pseudo)
RNA_pseudo <- FindVariableFeatures(RNA_pseudo, nfeatures = 10000)
RNA_pseudo[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo, pattern = "^MT-")
RNA_pseudo <- ScaleData(RNA_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_pseudo <- RunPCA(RNA_pseudo, npcs = 10)
RNA_pseudo <- RunUMAP(RNA_pseudo, dims = 1:10, reduction = 'pca',n.neighbors = 10L)
RNA_pseudo <- FindNeighbors(RNA_pseudo, dims = 1:10, reduction = "pca")
RNA_pseudo <- FindClusters(RNA_pseudo, resolution = 1)
pdf("./scRNA/RNA_pseudo_8cell.pdf")
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
dev.off()

label <- rownames(as.data.frame(RNA_pseudo@active.ident))
label <- as.data.frame(label)
label[,2] <- RNA_pseudo@meta.data$seurat_clusters

RNA_count <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_count <- RNA_count[which(rowSums(RNA_count) >0),]
toti <- read.csv("./totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("./pluripotency.markers.csv", sep = ",", header = F)
ZGA <- read.xlsx("/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/MajorZGA_genes.xlsx")
ZGA <- na.omit(ZGA$Gene)
RNA_ZGA <- RNA_count[which(rownames(RNA_count) %in% ZGA),]
RNA_toti <- RNA_count[which(rownames(RNA_count) %in% toti$V1),]
RNA_pluri <- RNA_count[which(rownames(RNA_count) %in% pluri$V1),]
RNA_ZGA <- colSums(RNA_ZGA)/620
RNA_toti <- colSums(RNA_toti)/215
RNA_pluri <- colSums(RNA_pluri)/38

label$ZGA <- RNA_ZGA
label$toti <- RNA_toti
label$pluri <- RNA_pluri
pdf("./scRNA/RNA_8cell_toti_pluri_average.pdf")
ggplot(label,mapping = aes(x = V2, y = pluri))+
  geom_violin(aes(fill =V2), trim = T) + 
  geom_boxplot(width = 0.05)+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0","#AC73EF"))+
  theme(legend.position = "none")+
  theme_bw()
dev.off()

RNA_pseudo@meta.data$zga.ident <- label$ZGA
RNA_pseudo@meta.data$zga.ident <- label$toti
RNA_pseudo@meta.data$zga.ident <- label$pluri

saveRDS(RNA_pseudo, "./scRNA/RNA_pseudo_8cell.rds")
############ hcluster 8cell stage ################
RNA_counts <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_counts <- as.matrix(RNA_counts)
RNA_counts <- t(RNA_counts)
out.dist <- dist(RNA_counts, method = "euclidean")
out.hclust <- hclust(out.dist, method = "ward.D2")
out.id <- cutree(out.hclust, k=2)
#identify late 2-cell as a cluster
pdf("./scRNA/RNA_8cell_pseudo.hcluster.pdf")
plot(out.hclust)
dev.off()

############################################# H3K4me3 (as an example) ####################################################
setwd("./h3k4me3")
h3k4me3 <- readRDS("./h3k4me3_allpeaks.rds")
h3k4me3 <- subset(h3k4me3,cells = colnames(h3k4me3[,grep("8cell", h3k4me3@meta.data$group.ident)]))
h3k4me3 <- RunTFIDF(h3k4me3)
h3k4me3 <- FindTopFeatures(h3k4me3, min.cutoff = 'q0')
h3k4me3 <- RunSVD(
  object = h3k4me3,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
h3k4me3 <- RunUMAP(object = h3k4me3, reduction = 'lsi', n.neighbors = 20L, dims = c(2:15))#或可调试reduction
h3k4me3 <- FindNeighbors(object = h3k4me3, reduction = 'lsi', dims = c(2:15))
h3k4me3 <- FindClusters(object = h3k4me3, algorithm = 3, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
h3k4me3_count <- as.matrix(h3k4me3@assays[["peaks"]]@counts)
depth <- apply(h3k4me3_count, 2, sum)
h3k4me3$depth <- log10(depth)

pdf("./h3k4me3_8cell.pdf")
DepthCor(h3k4me3)
DimPlot(object = h3k4me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(h3k4me3, reduction = "umap", group.by = "group.ident",label = TRUE, pt.size = 2)
FeaturePlot(h3k4me3, reduction = "umap", features = 'depth')
dev.off()

geneActivity <- read.csv("./h3k4me3_allstage_ActivitygeneActivityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(h3k4me3)]
h3k4me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k4me3) <- "ACTIVITY"
h3k4me3 <- FindVariableFeatures(h3k4me3, nfeatures = 10000)
h3k4me3 <- NormalizeData(h3k4me3)
h3k4me3 <- ScaleData(h3k4me3)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k4me3,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        #features = top500$gene,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, k.score = 10)
#Annotate scATAC-seq cells via label transfer

#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = h3k4me3[["lsi"]], dims = 2:ncol(h3k4me3[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k4me3_pre1 <- ID.predict
length(unique(h3k4me3_pre1$predicted.id))
#h3k4me3_pre1 <- h3k4me3_pre1[which(h3k4me3_pre1$prediction.score.max >0.3),]

write.table(h3k4me3_pre1, "../h3k4me3/h3k4me3_8cell_pseudo_inte.txt", quote = F, sep = " ")
saveRDS(h3k4me3, "../h3k4me3/h3k4me3_8cell_inte.rds")
saveRDS(transfer.anchors, "../h3k4me3/h3k4me3_8cell_transfer_anchors.rds")

####################################### only retain cells with prediction.score.max >0.2 across all modalities #######################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_8cell_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_8cell_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_8cell_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_8cell_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_8cell_h3k27ac_pseudo.txt", sep = " ")

h3k4me1_pre1 <- h3k4me1_pre1[which(h3k4me1_pre1$prediction.score.max >0.2),]
h3k4me3_pre1 <- h3k4me3_pre1[which(h3k4me3_pre1$prediction.score.max >0.2),]
h3k27ac_pre1 <- h3k27ac_pre1[which(h3k27ac_pre1$prediction.score.max >0.2),]
h3k36me3_pre1 <- h3k36me3_pre1[which(h3k36me3_pre1$prediction.score.max >0.2),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1 $predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1 [which(cotacit_pre1 $predicted.id %in% inte),]

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))
inte
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1$predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1$predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1$predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1$predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$predicted.id %in% inte),]

h3k4me1.sh <- paste0("cp ./bam/", h3k4me1_pre1$X, " ./8cell/", h3k4me1_pre1$predicted.id, "_inte.bam")
write.table(h3k4me1.sh, "./h3k4me1/h3k4me1_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k4me3.sh <- paste0("cp ./bam/", h3k4me3_pre1$X, " ./8cell/", h3k4me3_pre1$predicted.id, "_inte.bam")
write.table(h3k4me3.sh, "./h3k4me3/h3k4me3_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k27ac.sh <- paste0("cp ./bam/", h3k27ac_pre1$X, " ./8cell/", h3k27ac_pre1$predicted.id, "_inte.bam")
write.table(h3k27ac.sh, "./h3k27ac/h3k27ac_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k36me3.sh <- paste0("cp ./bam/", h3k36me3_pre1$X, " ./8cell/", h3k36me3_pre1$predicted.id, "_inte.bam")
write.table(h3k36me3.sh, "./h3k36me3/h3k36me3_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)


################################### annotation chromatin states #####################################3
######## 8cell_1 #########
C2_bed <- read.csv("./8cell_1/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 2),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 3),10] <- "Un"
C2_bed[which(C2_bed$V4 == 4),10] <- "Es"
C2_bed[which(C2_bed$V4 == 5),10] <- "Es"
C2_bed[which(C2_bed$V4 == 6),10] <- "Es"
C2_bed[which(C2_bed$V4 == 7),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 8),10] <- "Un"
C2_bed[which(C2_bed$V4 == 9),10] <- "K27"
C2_bed[which(C2_bed$V4 == 10),10] <- "K27"
C2_bed[which(C2_bed$V4 == 11),10] <- "He"
C2_bed[which(C2_bed$V4 == 12),10] <- "K9"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_1/8cell_1_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_2/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "He"
C2_bed[which(C2_bed$V4 == 2),10] <- "Un"
C2_bed[which(C2_bed$V4 == 3),10] <- "Un"
C2_bed[which(C2_bed$V4 == 4),10] <- "K27"
C2_bed[which(C2_bed$V4 == 5),10] <- "Un"
C2_bed[which(C2_bed$V4 == 6),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 7),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 8),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 9),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 10),10] <- "Es"
C2_bed[which(C2_bed$V4 == 11),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ps"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_2/8cell_2_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_3/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 2),10] <- "Un"
C2_bed[which(C2_bed$V4 == 3),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 4),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 5),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 6),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 7),10] <- "Es"
C2_bed[which(C2_bed$V4 == 8),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 9),10] <- "K27"
C2_bed[which(C2_bed$V4 == 10),10] <- "Un"
C2_bed[which(C2_bed$V4 == 11),10] <- "Un"
C2_bed[which(C2_bed$V4 == 12),10] <- "He"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_3/8cell_3_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_4/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "He"
C2_bed[which(C2_bed$V4 == 2),10] <- "Un"
C2_bed[which(C2_bed$V4 == 3),10] <- "K9"
C2_bed[which(C2_bed$V4 == 4),10] <- "K9"
C2_bed[which(C2_bed$V4 == 5),10] <- "K9"
C2_bed[which(C2_bed$V4 == 6),10] <- "Un"
C2_bed[which(C2_bed$V4 == 7),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 8),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 9),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 10),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 11),10] <- "Es"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ps"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_4/8cell_4_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_5/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "K27"
C2_bed[which(C2_bed$V4 == 2),10] <- "Un"
C2_bed[which(C2_bed$V4 == 3),10] <- "He"
C2_bed[which(C2_bed$V4 == 4),10] <- "K9"
C2_bed[which(C2_bed$V4 == 5),10] <- "Un"
C2_bed[which(C2_bed$V4 == 6),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 7),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 8),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 9),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 10),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 11),10] <- "Es"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ps"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_5/8cell_5_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_6/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "K27"
C2_bed[which(C2_bed$V4 == 2),10] <- "He"
C2_bed[which(C2_bed$V4 == 3),10] <- "Un"
C2_bed[which(C2_bed$V4 == 4),10] <- "K9"
C2_bed[which(C2_bed$V4 == 5),10] <- "He"
C2_bed[which(C2_bed$V4 == 6),10] <- "Un"
C2_bed[which(C2_bed$V4 == 7),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 8),10] <- "Un"
C2_bed[which(C2_bed$V4 == 9),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 10),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 11),10] <- "Es"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ps"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_6/8cell_6_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_7/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 2),10] <- "Es"
C2_bed[which(C2_bed$V4 == 3),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 4),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 5),10] <- "Es"
C2_bed[which(C2_bed$V4 == 6),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 7),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 8),10] <- "Un"
C2_bed[which(C2_bed$V4 == 9),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 10),10] <- "Un"
C2_bed[which(C2_bed$V4 == 11),10] <- "K9"
C2_bed[which(C2_bed$V4 == 12),10] <- "He"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_7/8cell_7_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")

C2_bed <- read.csv("./8cell_8/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Es"
C2_bed[which(C2_bed$V4 == 2),10] <- "Es"
C2_bed[which(C2_bed$V4 == 3),10] <- "Es"
C2_bed[which(C2_bed$V4 == 4),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 5),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 6),10] <- "Un"
C2_bed[which(C2_bed$V4 == 7),10] <- "Un"
C2_bed[which(C2_bed$V4 == 8),10] <- "K27"
C2_bed[which(C2_bed$V4 == 9),10] <- "He"
C2_bed[which(C2_bed$V4 == 10),10] <- "K9"
C2_bed[which(C2_bed$V4 == 11),10] <- "K9"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ps"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_8/8cell_8_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")










