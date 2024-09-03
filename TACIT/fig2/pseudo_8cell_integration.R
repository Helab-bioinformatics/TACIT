setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_all")

###########################################################################################################
############################### plote for 8cell_all chromHMM ##############################################
library(pheatmap)
C2 <- read.csv("./learnmodel_inte/emissions_12_chrhmm.txt", sep="\t")
rownames(C2) <- C2[,1]
C2 <- C2[,-1]
dat <- as.matrix(C2)
colnames(dat) <- colnames(C2)
rownames(dat) <- rownames(C2)
dat[dat<=0.01] <- 0
dat <- log10(dat+1)
pdf("8cell_intersect_all_nochrM_emissions_10_chrhmm.pdf")
p <- pheatmap(dat, 
              cluster_rows =F, cluster_cols = F, 
              show_rownames = T, scale ="none",
              show_colnames = T,
              clustering_method="ward.D2",
              clustering_distance_rows = "correlation",
              color = colorRampPalette(c("#D2EAF9","#277FC1","#2172B7"))(500))
#clustering_method should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
p
dev.off()

############################## annotation chromatin states ########################################
C2_bed <- read.csv("./learnmodel_inte/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 2),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 3),10] <- "Es"
C2_bed[which(C2_bed$V4 == 4),10] <- "Es"
C2_bed[which(C2_bed$V4 == 5),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 6),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 7),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 8),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 9),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 10),10] <- "Un"
C2_bed[which(C2_bed$V4 == 11),10] <- "K27"
C2_bed[which(C2_bed$V4 == 12),10] <- "He"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "./8cell_all_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")


C2_bed_E <- C2_bed[grep("Ep|Es|Ew", C2_bed$V10),]
write.table(C2_bed_E, "./8cell_all_chrhmm_enhancer.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_P <- C2_bed[grep("Pw|Ps", C2_bed$V10),]
write.table(C2_bed_P, "./8cell_all_chrhmm_promoter.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_T <- C2_bed[grep("Tw|Ts", C2_bed$V10),]
write.table(C2_bed_T, "./8cell_all_chrhmm_transcription.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_M <- C2_bed[grep("M", C2_bed$V10),]
write.table(C2_bed_M, "./8cell_all_chrhmm_multivalent.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K9 <- C2_bed[grep("K9", C2_bed$V10),]
write.table(C2_bed_K9, "./8cell_all_chrhmm_h3k9me3.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("K27", C2_bed$V10),]
write.table(C2_bed_K27, "./8cell_all_chrhmm_polycomb.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("He", C2_bed$V10),]
write.table(C2_bed_K27, "./8cell_all_chrhmm_heterochromtin.bed", quote = F, col.names = F, row.names = F, sep = "\t")

####### change colors for bed file ##############
C1_bed <- read.csv("./8cell_all_mm10_chrhmm_dense_anno.bed", sep = "\t", header = F)
C1_bed[grep("Un", C1_bed$V4),9] <- "209,211,212"
C1_bed[grep("Mu", C1_bed$V4),9] <- "170,104,55"
C1_bed[grep("Pw", C1_bed$V4),9] <- "117,204,158"
C1_bed[grep("Ps", C1_bed$V4),9] <- "11,158,91"
C1_bed[grep("Ew", C1_bed$V4),9] <- "249,209,160"
C1_bed[grep("Es", C1_bed$V4),9] <- "249,156,8"
C1_bed[grep("Tw", C1_bed$V4),9] <- "214,107,123"
C1_bed[grep("Ts", C1_bed$V4),9] <- "198,32,48"
C1_bed[grep("K27", C1_bed$V4),9] <- "146,234,244"
C1_bed[grep("K9", C1_bed$V4),9] <- "21,162,206"
C1_bed[grep("He", C1_bed$V4),9] <- "109,110,113"
write.table(C1_bed, "./8cell_all_mm10_chrhmm_dense_anno_color.bed", quote = F, col.names = F, row.names = F, sep = "\t")


#######################################################################################################################
####################################### pseudo-cell integration #######################################################
library(cisTopic)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(dplyr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(edgeR)
library(openxlsx)
set.seed(1)

############################################# RNA ########################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/")

#######################################################################################################################
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
        panel.grid.minor = element_blank())
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

############################################# H3K4me3 ####################################################
setwd("./h3k4me3")
h3k4me3 <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allpeaks.rds")
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

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allstage_ActivitygeneActivityrix.txt", sep = " ",check.names=F)
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

######################################### H3K4me1 #############################################
############H3K4me1的4细胞数目少且reads数多，对应效果好，不需要构建pseudo-cells
setwd("../h3k4me1/")
set.seed(1)
h3k4me1_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/h3k4me1.5kb.countMatrix.HX.txt", sep ="\t") 
rownames(h3k4me1_counts) <- h3k4me1_counts[,1]
h3k4me1_counts <- h3k4me1_counts[,-1]
h3k4me1_counts <- h3k4me1_counts[,grep("8cell",colnames(h3k4me1_counts))] 
h3k4me1 <- CreateSeuratObject(counts = h3k4me1_counts,assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = 'q0')
h3k4me1 <- RunSVD(object = h3k4me1, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

DefaultAssay(h3k4me1) <- "peaks"
pdf("h3k4me1_8cell_5kbbin.pdf")
DepthCor(h3k4me1)
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:10), n.neighbors = 10L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) 
label <- rownames(as.data.frame(h3k4me1@active.ident))
label <- as.data.frame(label)
label[grep("8cell",label[,1]),2] <- "8cell"
h3k4me1@meta.data$group.ident <- as.factor(label[,2])
DimPlot(h3k4me1, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off()

#geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me1/h3k4me1_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me1/h3k4me1_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(h3k4me1)]
h3k4me1[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k4me1) <- "ACTIVITY"
h3k4me1 <- FindVariableFeatures(h3k4me1, nfeatures = 10000)
h3k4me1 <- NormalizeData(h3k4me1)
h3k4me1 <- ScaleData(h3k4me1)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k4me1,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)

#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = h3k4me1[["lsi"]], dims = 2:ncol(h3k4me1[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k4me1_pre1 <- ID.predict
length(unique(h3k4me1_pre1$predicted.id))
#h3k4me1_pre1 <- h3k4me1_pre1[which(h3k4me1_pre1$prediction.score.max >0.2),]
write.table(h3k4me1_pre1, "../h3k4me1/h3k4me1_8cell_pseudo.txt", quote = F, sep = " ")

saveRDS(h3k4me1, "../h3k4me1/h3k4me1_8cell_integration.rds")
saveRDS(transfer.anchors, "../h3k4me1/h3k4me1_8cell_transfer_anchors.rds")

######################################### H3K27ac #######################################
setwd("../h3k27ac/")
set.seed(2)
h3k27ac_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/h3k27ac.5kb.best.countMatrix.txt", sep =" ") 
h3k27ac_counts <- h3k27ac_counts[,grep("8cell",colnames(h3k27ac_counts))] 
h3k27ac <- CreateSeuratObject(counts = h3k27ac_counts,assay = 'peaks',project = 'h3k27ac',min.cells = 1)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = 'q0')
h3k27ac <- RunSVD(object = h3k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k27ac/h3k27ac_allStage_peaks_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,grep("8cell", colnames(geneActivity))]

dat <- data.frame(col=colnames(geneActivity), min="min")
dat_1 <- dat[1:41,]
dat_2 <- dat[42:89,]
dat_1 <- separate(dat_1, col, c("v1","v2", "v3"), sep = "-")
dat_1 <- unite(dat_1,col, c("v1","v2", "v3"), sep = ".")
dat <- rbind(dat_1, dat_2)
colnames(geneActivity) <- dat$col
geneActivity <- geneActivity[,colnames(h3k27ac)]

h3k27ac[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k27ac) <- "ACTIVITY"
h3k27ac <- FindVariableFeatures(h3k27ac, nfeatures = 10000)
h3k27ac <- NormalizeData(h3k27ac)
h3k27ac <- ScaleData(h3k27ac)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k27ac,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)

#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = h3k27ac[["lsi"]], dims = 2:ncol(h3k27ac[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k27ac_pre1 <- ID.predict
length(unique(h3k27ac_pre1$predicted.id))
#h3k27ac_pre1 <- h3k27ac_pre1[which(h3k27ac_pre1$prediction.score.max >0.3),]
write.table(h3k27ac_pre1, "./h3k27ac_8cell_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k27ac, "h3k27ac_8cell_integration.rds")
saveRDS(transfer.anchors, "h3k27ac_8cell_transfer_anchors.rds")


############################################### H3K36me3 #########################################
#rm(list = ls())
setwd("../h3k36me3/")
set.seed(3)
h3k36me3 <- readRDS("./h3k36me3_5kb.rds")
h3k36me3 <- subset(h3k36me3,cells = colnames(h3k36me3[,grep("8cell", h3k36me3@meta.data$group.ident)]))
pdf("h3k36me3_8cell_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 
####目前的UMAP不合适，有batch effect，会影响trajectory，所以重新降维
h3k36me3 <- RunUMAP(object = h3k36me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
h3k36me3 <- FindNeighbors(object = h3k36me3,reduction = 'lsi',dims = c(2:15))
h3k36me3 <- FindClusters(object = h3k36me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
pdf("h3k36me3_8cell_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k36me3/h3k36me3_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,grep("8cell", colnames(geneActivity))]
geneActivity <- geneActivity[,colnames(h3k36me3)]
h3k36me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k36me3) <- "ACTIVITY"
h3k36me3 <- FindVariableFeatures(h3k36me3, nfeatures = 10000)
h3k36me3 <- NormalizeData(h3k36me3)
h3k36me3 <- ScaleData(h3k36me3)

#RNA <- readRDS("../RNA_8cell.rds")
#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k36me3,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = h3k36me3[["lsi"]], dims = 2:ncol(h3k36me3[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k36me3_pre1 <- ID.predict
length(unique(h3k36me3_pre1$predicted.id))
#h3k36me3_pre1 <- h3k36me3_pre1[which(h3k36me3_pre1$prediction.score.max >0.2),]
write.table(h3k36me3_pre1, "./h3k36me3_8cell_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k36me3, "h3k36me3_8cell_integration.rds")
saveRDS(transfer.anchors, "h3k36me3_8cell_transfer_anchors.rds")

########################################### H3K27me3 #############################################
##借助coTACIT数据
setwd("../cotacit_data/")
cotacit <- readRDS("./cotacit_h3k27ac.rds")
cotacit_8cell <- subset(cotacit,cells = colnames(cotacit[,grep("8cell", cotacit@meta.data$group.ident)]))
geneActivity <- read.csv("./cotacit_h3k27ac_allstage_tss5kb_ActivityMatrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(cotacit_8cell)]
cotacit_8cell[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(cotacit_8cell) <- "ACTIVITY"
cotacit_8cell <- FindVariableFeatures(cotacit_8cell, nfeatures = 10000)
cotacit_8cell <- NormalizeData(cotacit_8cell)
cotacit_8cell <- ScaleData(cotacit_8cell)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = cotacit_8cell,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = cotacit_8cell[["lsi"]], dims = 2:ncol(cotacit_8cell[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_pre1 <- ID.predict
length(unique(cotacit_pre1$predicted.id))
#cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]
write.table(cotacit_pre1, "./cotacit_8cell_h3k27ac_pseudo.txt", quote = F, sep = " ")
saveRDS(cotacit_8cell, "cotacit_8cell_h3k27ac_integration.rds")
saveRDS(transfer.anchors, "cotacit_8cell_h3k27ac_transfer_anchors.rds")


####################################### all #######################################################
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

write.table(h3k4me1_pre1, "./h3k4me1_8cell_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k4me3_pre1, "./h3k4me3_8cell_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k27ac_pre1, "./h3k27ac_8cell_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k36me3_pre1, "./h3k36me3_8cell_pseudo_select.txt", quote = F, sep = " ")
write.table(cotacit_pre1, "./cotacit_8cell_pseudo_select.txt", quote = F, sep = " ")
##由于想尽可能多地保留整合后的细胞，所以整合到同一个RNA细胞的多个细胞也分开并保留，因此在predicted.id后面还增加了一个编号，这一步在excel中操作(按照predicted.id和score值进行sort后编号)。
h3k4me1_pre1 <- read.csv("./h3k4me1_8cell_pseudo_select.txt", sep = "\t")
h3k4me3_pre1 <- read.csv("./h3k4me3_8cell_pseudo_select.txt", sep = "\t")
h3k27ac_pre1 <- read.csv("./h3k27ac_8cell_pseudo_select.txt", sep = "\t")
h3k36me3_pre1 <- read.csv("./h3k36me3_8cell_pseudo_select.txt", sep = "\t")
cotacit_pre1 <- read.csv("./cotacit_8cell_pseudo_select.txt", sep = "\t")

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))
inte
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1$predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1$predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1$predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1$predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$predicted.id %in% inte),]

h3k4me1.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/", h3k4me1_pre1$X, " ./8cell/", h3k4me1_pre1$predicted.id, "_inte.bam")
write.table(h3k4me1.sh, "./h3k4me1/h3k4me1_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", h3k4me3_pre1$X, " ./8cell/", h3k4me3_pre1$predicted.id, "_inte.bam")
write.table(h3k4me3.sh, "./h3k4me3/h3k4me3_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k27ac_pre1_1 <- h3k27ac_pre1[grep("2238",h3k27ac_pre1$X),]
h3k27ac_pre1_2 <- h3k27ac_pre1[grep("2236",h3k27ac_pre1$X),]
h3k27ac_pre1_1 <- separate(h3k27ac_pre1_1, X, c("x1", "x2", "x3", "x4", "x5"), sep = "\\.")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, id, c("x1", "x2", "x3"), sep = "-")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, X, c("id", "x4", "x5"), sep = ".")
h3k27ac_pre1 <- rbind(h3k27ac_pre1_1, h3k27ac_pre1_2)
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", h3k27ac_pre1$X, " ./8cell/", h3k27ac_pre1$predicted.id, "_inte.bam")
write.table(h3k27ac.sh, "./h3k27ac/h3k27ac_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", h3k36me3_pre1$X, " ./8cell/", h3k36me3_pre1$predicted.id, "_inte.bam")
write.table(h3k36me3.sh, "./h3k36me3/h3k36me3_8cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

cotacit_pre1$ID <- cotacit_pre1$X
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_combined_all_h3k9me3.mm10_rmdup.bam")
cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2, " ./8cell/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k27me3/cotacit_8cell_linked_h3k27me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3, " ./8cell/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k9me3/cotacit_8cell_linked_h3k9me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(cotacit_pre1[,c("predicted.id","prediction.score.max", "ID2", "ID3")], "./cotacit_data/integration/cotacit_8cell_integration.txt")


############################################# cotacit + tacit h3k37me3 #########################################

############ cotacit h3k27me3 ################
library(Signac)
h3k27me3 <- readRDS("./h3k27me3_all stage_peaks.rds")
cotacit_h3k27me3 <- readRDS("./cotacit_data/cotacit_h3k27me3_peaks.rds")

h3k27me3$dataset <- "tacit"
cotacit_h3k27me3$dataset <- "cotacit"

h3k27me3_8cell <- subset(h3k27me3,cells = colnames(h3k27me3)[grep("8cell", colnames(h3k27me3))])
cotacit_h3k27me3_8cell <- subset(cotacit_h3k27me3,cells = colnames(cotacit_h3k27me3)[grep("8cell", colnames(cotacit_h3k27me3))])
merge <- merge(h3k27me3_8cell, cotacit_h3k27me3_8cell)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_8cell.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k27me3_8cell, cotacit_h3k27me3_8cell),
  anchor.features = rownames(h3k27me3_8cell),
  reduction = "rlsi",
  dims = 2:30)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merge[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 2:30,
  k.weight=20)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_integration_8cell.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k27me3_8cell,
  query = cotacit_h3k27me3_8cell,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k27me3_8cell@meta.data$ID <- colnames(h3k27me3_8cell)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k27me3_8cell$ID,
                                   weight.reduction = cotacit_h3k27me3_8cell[["lsi"]], dims = 2:ncol(cotacit_h3k27me3_8cell[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k27me3_pre1 <- ID.predict
length(unique(cotacit_h3k27me3_pre1$predicted.id))

write.table(cotacit_h3k27me3_pre1, "./cotacit_data/h3k27me3_8cell_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k27me3_8cell_cotacit_tacit_inte_transfer_anchors.rds")

############################################# cotacit + tacit h3k9me3 #########################################

############ cotacit h3k9me3 ################
cotacit_h3k9me3 <- readRDS("./cotacit_data/cotacit_h3k9me3_10kbbin.rds")
h3k9me3 <- readRDS("./cotacit_data/h3k9me3_10kbbin.rds")
h3k9me3$dataset <- "tacit"
cotacit_h3k9me3$dataset <- "cotacit"

h3k9me3_8cell <- subset(h3k9me3,cells = colnames(h3k9me3)[grep("8cell", colnames(h3k9me3))])
cotacit_h3k9me3_8cell <- subset(cotacit_h3k9me3,cells = colnames(cotacit_h3k9me3)[grep("8cell", colnames(cotacit_h3k9me3))])
merge <- merge(h3k9me3_8cell, cotacit_h3k9me3_8cell)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_8cell.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k9me3_8cell, cotacit_h3k9me3_8cell),
  anchor.features = rownames(h3k9me3_8cell),
  reduction = "rlsi",
  dims = 2:30)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merge[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  k.weight=15)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_integration_8cell.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k9me3_8cell,
  query = cotacit_h3k9me3_8cell,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k9me3_8cell@meta.data$ID <- colnames(h3k9me3_8cell)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k9me3_8cell$ID,
                                   weight.reduction = cotacit_h3k9me3_8cell[["lsi"]], dims = 2:ncol(cotacit_h3k9me3_8cell[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k9me3_pre1 <- ID.predict
length(unique(cotacit_h3k9me3_pre1$predicted.id))

write.table(cotacit_h3k9me3_pre1, "./cotacit_data/h3k9me3_8cell_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k9me3_8cell_cotacit_tacit_inte_transfer_anchors.rds")

cotacit_pre1$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1$ID3, rownames(cotacit_h3k9me3_pre1))]
k27me3.sh <- paste0("samtools merge ./h3k27me3_tacit+cotacit/8cell/", cotacit_pre1$predicted.id, "_inte.bam ", 
                    "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2,
                    " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1$k27_tacit)
k9me3.sh <- paste0("samtools merge ./h3k9me3_tacit+cotacit/8cell/", cotacit_pre1$predicted.id, "_inte.bam ", 
                   "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3,
                   " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1$k9_tacit)
write.table(k27me3.sh, "./cotacit_data/integration/h3k27me3_8cell_merge.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./cotacit_data/integration/h3k9me3_8cell_merge.sh", sep = "\t", quote = F, col.names = F, row.names = F)

write.table(cotacit_pre1, "./cotacit_data/cotacit_tacit_inte_8cell.txt", quote = F, sep = " ")


################################## 为了用尽可能多的细胞构建chromHMM，基本保留所有整合的TACIT和cotacit细胞 ###########################################
###### 这里和pseudo_8cell_meta.R做的事情类似 #####
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

#### h3k4me1 ####
dir.sh <- paste0("mkdir ", inte)
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1 $predicted.id %in% inte),]
sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/", rownames(h3k4me1_pre1), 
             " ./", h3k4me1_pre1$predicted.id, "/")
h3k4me1.sh <- c(dir.sh, sh)
write.table(h3k4me1.sh, "./04_chromHMM_8cell/h3k4me1_8cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k4me3 ####
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", rownames(h3k4me3_pre1), 
                     " ./", h3k4me3_pre1$predicted.id, "/")
write.table(h3k4me3.sh, "./04_chromHMM_8cell/h3k4me3_8cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k27ac ####
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k27ac_pre1$X <- rownames(h3k27ac_pre1)
h3k27ac_pre1_1 <- h3k27ac_pre1[grep("2238",h3k27ac_pre1$X),]
h3k27ac_pre1_2 <- h3k27ac_pre1[grep("2236",h3k27ac_pre1$X),]
h3k27ac_pre1_1 <- separate(h3k27ac_pre1_1, X, c("x1", "x2", "x3", "x4", "x5"), sep = "\\.")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, id, c("x1", "x2", "x3"), sep = "-")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, X, c("id", "x4", "x5"), sep = ".")
h3k27ac_pre1 <- rbind(h3k27ac_pre1_1, h3k27ac_pre1_2)
rownames(h3k27ac_pre1) <- h3k27ac_pre1$X
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", rownames(h3k27ac_pre1), 
                     " ./", h3k27ac_pre1$predicted.id, "/")
write.table(h3k27ac.sh, "./04_chromHMM_8cell/h3k27ac_8cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k36me3 ####
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", rownames(h3k36me3_pre1), 
                      " ./", h3k36me3_pre1$predicted.id, "/")
write.table(h3k36me3.sh, "./04_chromHMM_8cell/h3k36me3_8cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)


#### cotacit ####
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_8cell_h3k27ac_pseudo.txt", sep = " ")
cotacit_h3k27me3_pre1 <- read.csv("./cotacit_data/h3k27me3_8cell_cotacit_tacit_inte.txt", sep = " ")
cotacit_h3k9me3_pre1 <- read.csv("./cotacit_data/h3k9me3_8cell_cotacit_tacit_inte.txt", sep = " ")

cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1 $predicted.id %in% inte),]
cotacit_pre1$ID <- rownames(cotacit_pre1)
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_combined_all_h3k9me3.mm10_rmdup.bam") 
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
write.table(k27me3.sh, "./04_chromHMM_8cell/h3k27me3_8cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./04_chromHMM_8cell/h3k9me3_8cell_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)



##################################### plot for chromHMM #############################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/")

library(pheatmap)
C2 <- read.csv("./8cell_1/learnmodel/emissions_12_chrhmm.txt", sep="\t")
rownames(C2) <- C2[,1]
C2 <- C2[,-1]
dat <- as.matrix(C2)
colnames(dat) <- colnames(C2)
rownames(dat) <- rownames(C2)
dat[dat<=0.01] <- 0
dat <- log10(dat+1)
pdf("8cell_1_intersect_nochrM_emissions_12_chrhmm.pdf")
p <- pheatmap(dat, 
              cluster_rows =F, cluster_cols = F, 
              show_rownames = T, scale ="none",
              show_colnames = T,
              clustering_method="ward.D2",
              clustering_distance_rows = "correlation",
              color = colorRampPalette(c("#D2EAF9","#277FC1","#2172B7"))(500))
#clustering_method should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
p
dev.off()


################################### annotation chromatin states #####################################3
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/")
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










