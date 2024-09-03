

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

######################## 直接看single cell和pseudo-cells，能否分成ICM和TE potential，还是单细胞效果好，故ICM和TE分别合成pseudo-cells #########################
RNA_pseudo <- readRDS("./RNA_pseudo.rds")
RNA_pseudo <- subset(RNA_pseudo,cells = colnames(RNA_pseudo[,grep("blasto", colnames(RNA_pseudo))]))
RNA_pseudo <- NormalizeData(RNA_pseudo)
RNA_pseudo <- FindVariableFeatures(RNA_pseudo, nfeatures = 2000)
RNA_pseudo[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo, pattern = "^MT-")
RNA_pseudo <- ScaleData(RNA_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_pseudo <- RunPCA(RNA_pseudo, npcs = 30)
RNA_pseudo <- RunUMAP(RNA_pseudo, dims = 1:30, reduction = 'pca',n.neighbors = 10L)
RNA_pseudo <- FindNeighbors(RNA_pseudo, dims = 1:30, reduction = "pca")
RNA_pseudo <- FindClusters(RNA_pseudo, resolution = 1)

pdf("./scRNA/RNA_pseudo_blasto_markers.pdf")
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, group.by = "group.ident" ,label = TRUE, pt.size = 2)
FeaturePlot(object = RNA_pseudo,features = c("Sox2", "Nanog", "Pou5f1", "Klf4", "Otx2", "Gata4"),pt.size = 1.5,max.cutoff = 'q95', ncol = 2)
FeaturePlot(object = RNA_pseudo,features = c("Cdx2", "Eomes", "Gata2", "Id2", "Tbx2", "Elf5"),pt.size = 1.5,max.cutoff = 'q95', ncol = 2)
dev.off()


RNA <- readRDS("./scRNA/all_RNA_no_outlier.rds")
RNA <- subset(RNA,cells = colnames(RNA[,grep("blasto", colnames(RNA))]))
RNA <- NormalizeData(RNA)
RNA <- FindVariableFeatures(RNA, nfeatures = 3000)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = "^MT-")
RNA <- ScaleData(RNA, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA <- RunPCA(RNA, npcs = 50)
RNA <- RunUMAP(RNA, dims = 1:50, reduction = 'pca',n.neighbors = 5L)
RNA <- FindNeighbors(RNA, dims = 1:50, reduction = "pca")
RNA <- FindClusters(RNA, resolution = 0.5)

pdf("./scRNA/RNA_blasto_markers.pdf", width = 12)
DimPlot(RNA, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA, group.by = "group.ident" ,label = TRUE, pt.size = 2)
FeaturePlot(object = RNA,features = c("Nanog", "Klf4",  "Etv5", "Fgf4",  "Utf1", "Esrrb","Pdgfra"),pt.size = 1.5,max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = RNA,features = c("Cdx2", "Eomes", "Gata3", "Id2",  "Krt8", "Dppa1", "Krt18", "Fgfr2"),pt.size = 1.5,max.cutoff = 'q95', ncol = 3)
dev.off()

################################ ICM和TE分别合成pseudo-cells #########################
RNA_ICM <- subset(RNA,cells = colnames(RNA[,grep("0", RNA@meta.data[["seurat_clusters"]])]))
RNA_TE <- subset(RNA,cells = colnames(RNA[,grep("1", RNA@meta.data[["seurat_clusters"]])]))

####### ICM ##########
RNA_ICM <- NormalizeData(RNA_ICM)
RNA_ICM <- FindVariableFeatures(RNA_ICM, nfeatures = 3000)
RNA_ICM[["percent.mt"]] <- PercentageFeatureSet(RNA_ICM, pattern = "^MT-")
RNA_ICM <- ScaleData(RNA_ICM, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_ICM <- RunPCA(RNA_ICM, npcs = 10)
RNA_ICM <- RunUMAP(RNA_ICM, dims = 1:10, reduction = 'pca',n.neighbors = 5L)
RNA_ICM <- FindNeighbors(RNA_ICM, dims = 1:10, reduction = "pca")
RNA_ICM <- FindClusters(RNA_ICM, resolution = 0.5)

library(monocle3)
set.seed(1234)
data <- GetAssayData(RNA_ICM, assay = 'RNA', slot = 'counts')
label <- as.data.frame(colnames(RNA_ICM))
label[grep("earlyblast",label[,1]),2] <- "earlyblast"
label[grep("lateblast",label[,1]),2] <- "lateblast"
colnames(label) <- c("ID", "stage")
rownames(label) <- label$ID
cds <- new_cell_data_set(data,cell_metadata = label)
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(RNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('int.umap')
pdf("./scRNA/RNA_ICM_trajectory.pdf")
p1
p2
## 识别轨迹
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "stage", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE,  graph_label_size=3, cell_size = 1.2)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE,  graph_label_size=3, cell_size = 1.2)
dev.off()
pseudo_order <- sort(cds@principal_graph_aux@listData$UMAP$pseudotime)
pseudo_order <- as.data.frame(pseudo_order)
saveRDS(cds, "./scRNA/RNA_ICM_trajectory.rds")

RNA_counts <- as.matrix(RNA_ICM@assays[["RNA"]]@counts)
RNA_counts <- RNA_counts[,rownames(pseudo_order)]

#pesudo-cells, merge every 5 cells, for each stage one by one
dim(RNA_counts)#16619   118
dat2 <- matrix(1:382237,ncol=23)
id <- seq(from=1,to=115,by=5)
n=1
for (i in id) {
  m=i+4
  dat2[,n] <- rowSums(RNA_counts[,i:m])
  n=n+1
}
dat2 <- as.data.frame(dat2)
rownames(dat2) <- rownames(RNA_counts)
colnames(dat2) <- paste0("ICM_pseudo_", 1:ncol(dat2))
# 
# RNA_pseudo_ICM <- CreateSeuratObject(counts = dat3, project = "RNA", min.cells = 3, min.features = 3000)
# RNA_pseudo_ICM[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo_ICM, pattern = "^MT-")
# pdf("./scRNA/RNA_pseudo_ICM.pdf")
# VlnPlot(RNA_pseudo_ICM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
# plot1 <- FeatureScatter(RNA_pseudo_ICM, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(RNA_pseudo_ICM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 
# plot2
# RNA_pseudo_ICM <- NormalizeData(RNA_pseudo_ICM)
# RNA_pseudo_ICM <- FindVariableFeatures(RNA_pseudo_ICM, nfeatures = 8000)
# RNA_pseudo_ICM <- ScaleData(RNA_pseudo_ICM, vars.to.regress = c("nFeature_RNA", "percent.mito"))
# RNA_pseudo_ICM <- RunPCA(RNA_pseudo_ICM, npcs = 15)
# RNA_pseudo_ICM <- RunUMAP(RNA_pseudo_ICM, dims = 1:15, reduction = 'pca',n.neighbors = 5L)
# RNA_pseudo_ICM <- FindNeighbors(RNA_pseudo_ICM, dims = 1:15, reduction = "pca")
# RNA_pseudo_ICM <- FindClusters(RNA_pseudo_ICM, resolution = 0.5)
# DimPlot(RNA_pseudo_ICM, reduction = "umap", label = TRUE, pt.size = 2)
# DimPlot(RNA_pseudo_ICM, reduction = "pca", label = TRUE, pt.size = 2)
# dev.off()
# saveRDS(RNA_pseudo_ICM, "./scRNA/RNA_pseudo_ICM.rds")


####### TE ##########
RNA_TE <- NormalizeData(RNA_TE)
RNA_TE <- FindVariableFeatures(RNA_TE, nfeatures = 3000)
RNA_TE[["percent.mt"]] <- PercentageFeatureSet(RNA_TE, pattern = "^MT-")
RNA_TE <- ScaleData(RNA_TE, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_TE <- RunPCA(RNA_TE, npcs = 10)
RNA_TE <- RunUMAP(RNA_TE, dims = 1:10, reduction = 'pca',n.neighbors = 5L)
RNA_TE <- FindNeighbors(RNA_TE, dims = 1:10, reduction = "pca")
RNA_TE <- FindClusters(RNA_TE, resolution = 0.5)

library(monocle3)
set.seed(1234)
data <- GetAssayData(RNA_TE, assay = 'RNA', slot = 'counts')
label <- as.data.frame(colnames(RNA_TE))
label[grep("earlyblast",label[,1]),2] <- "earlyblast"
label[grep("lateblast",label[,1]),2] <- "lateblast"
colnames(label) <- c("ID", "stage")
rownames(label) <- label$ID
cds <- new_cell_data_set(data,cell_metadata = label)
cds <- preprocess_cds(cds, num_dim = 50)
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(RNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('int.umap')
pdf("./scRNA/RNA_TE_trajectory.pdf")
p1
p2
## 识别轨迹
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "stage", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE,  graph_label_size=3, cell_size = 1.2)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE,  graph_label_size=3, cell_size = 1.2)
dev.off()
pseudo_order <- sort(cds@principal_graph_aux@listData$UMAP$pseudotime)
pseudo_order <- as.data.frame(pseudo_order)
saveRDS(cds, "./scRNA/RNA_TE_trajectory.rds")

RNA_counts <- as.matrix(RNA_TE@assays[["RNA"]]@counts)
RNA_counts <- RNA_counts[,rownames(pseudo_order)]

#pesudo-cells, merge every 5 cells, for each stage one by one
dim(RNA_counts)#16619   97
dat3 <- matrix(1:315761,ncol=19)
id <- seq(from=1,to=95,by=5)
n=1
for (i in id) {
  m=i+4
  dat3[,n] <- rowSums(RNA_counts[,i:m])
  n=n+1
}
dat3 <- as.data.frame(dat3)
rownames(dat3) <- rownames(RNA_counts)
colnames(dat3) <- paste0("TE_pseudo_", 1:ncol(dat3))

# RNA_pseudo_TE <- CreateSeuratObject(counts = dat3, project = "RNA", min.cells = 3, min.features = 3000)
# RNA_pseudo_TE[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo_TE, pattern = "^MT-")
# pdf("./scRNA/RNA_pseudo_TE.pdf")
# VlnPlot(RNA_pseudo_TE, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
# plot1 <- FeatureScatter(RNA_pseudo_TE, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(RNA_pseudo_TE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 
# plot2
# RNA_pseudo_TE <- NormalizeData(RNA_pseudo_TE)
# RNA_pseudo_TE <- FindVariableFeatures(RNA_pseudo_TE, nfeatures = 8000)
# RNA_pseudo_TE <- ScaleData(RNA_pseudo_TE, vars.to.regress = c("nFeature_RNA", "percent.mito"))
# RNA_pseudo_TE <- RunPCA(RNA_pseudo_TE, npcs = 15)
# RNA_pseudo_TE <- RunUMAP(RNA_pseudo_TE, dims = 1:15, reduction = 'pca',n.neighbors = 5L)
# RNA_pseudo_TE <- FindNeighbors(RNA_pseudo_TE, dims = 1:15, reduction = "pca")
# RNA_pseudo_TE <- FindClusters(RNA_pseudo_TE, resolution = 0.5)
# DimPlot(RNA_pseudo_TE, reduction = "umap", label = TRUE, pt.size = 2)
# DimPlot(RNA_pseudo_TE, reduction = "pca", label = TRUE, pt.size = 2)
# dev.off()
# saveRDS(RNA_pseudo_TE, "./scRNA/RNA_pseudo_TE.rds")

############ combine ################
all <- cbind(dat2, dat3)
blasto_pseudo <- CreateSeuratObject(counts = all, project = "RNA", min.cells = 3, min.features = 3000)
blasto_pseudo[["percent.mt"]] <- PercentageFeatureSet(blasto_pseudo, pattern = "^MT-")
pdf("./scRNA/blasto_pseudo.pdf", width = 12)
VlnPlot(blasto_pseudo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
plot1 <- FeatureScatter(blasto_pseudo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(blasto_pseudo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
blasto_pseudo <- NormalizeData(blasto_pseudo)
blasto_pseudo <- FindVariableFeatures(blasto_pseudo, nfeatures = 3000)
blasto_pseudo <- ScaleData(blasto_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
blasto_pseudo <- RunPCA(blasto_pseudo, npcs = 30)
blasto_pseudo <- RunUMAP(blasto_pseudo, dims = 1:30, reduction = 'pca',n.neighbors = 5L)
blasto_pseudo <- FindNeighbors(blasto_pseudo, dims = 1:30, reduction = "pca")
blasto_pseudo <- FindClusters(blasto_pseudo, resolution = 0.5)
DimPlot(blasto_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
FeaturePlot(object = blasto_pseudo,features = c("Nanog", "Klf4",  "Etv5", "Fgf4",  "Utf1", "Esrrb","Pdgfra"),pt.size = 1.5,max.cutoff = 'q95', ncol = 3)
FeaturePlot(object = blasto_pseudo,features = c("Cdx2", "Eomes", "Gata3", "Id2",  "Krt8", "Dppa1", "Krt18", "Fgfr2"),pt.size = 1.5,max.cutoff = 'q95', ncol = 3)
dev.off()
saveRDS(blasto_pseudo, "./scRNA/blasto_pseudo.rds")

############################################# H3K4me3 ####################################################
setwd("./h3k4me3")
h3k4me3 <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allpeaks.rds")
h3k4me3 <- subset(h3k4me3,cells = colnames(h3k4me3[,grep("blasto", h3k4me3@meta.data$group.ident)]))
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

pdf("./h3k4me3_blasto.pdf")
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
transfer.anchors <- FindTransferAnchors(reference = blasto_pseudo,
                                        query = h3k4me3,
                                        features = VariableFeatures(object = blasto_pseudo),
                                        #features = top500$gene,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, k.score = 10)
#Annotate scATAC-seq cells via label transfer

#transfer ID
blasto_pseudo@meta.data$ID <- colnames(blasto_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = blasto_pseudo$ID,
                                   weight.reduction = h3k4me3[["lsi"]], dims = 2:ncol(h3k4me3[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k4me3_pre1 <- ID.predict
length(unique(h3k4me3_pre1$predicted.id))
#h3k4me3_pre1 <- h3k4me3_pre1[which(h3k4me3_pre1$prediction.score.max >0.3),]

write.table(h3k4me3_pre1, "../h3k4me3/h3k4me3_blasto_pseudo_inte.txt", quote = F, sep = " ")
saveRDS(h3k4me3, "../h3k4me3/h3k4me3_blasto_inte.rds")
saveRDS(transfer.anchors, "../h3k4me3/h3k4me3_blasto_transfer_anchors.rds")

######################################### H3K4me1 #############################################
############H3K4me1的4细胞数目少且reads数多，对应效果好，不需要构建pseudo-cells
setwd("../h3k4me1/")
set.seed(1)
h3k4me1_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/h3k4me1.5kb.countMatrix.HX.txt", sep ="\t") 
rownames(h3k4me1_counts) <- h3k4me1_counts[,1]
h3k4me1_counts <- h3k4me1_counts[,-1]
h3k4me1_counts <- h3k4me1_counts[,grep("blasto",colnames(h3k4me1_counts))] 
h3k4me1 <- CreateSeuratObject(counts = h3k4me1_counts,assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = 'q0')
h3k4me1 <- RunSVD(object = h3k4me1, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

DefaultAssay(h3k4me1) <- "peaks"
pdf("h3k4me1_blasto_5kbbin.pdf")
DepthCor(h3k4me1)
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:10), n.neighbors = 10L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) 
label <- rownames(as.data.frame(h3k4me1@active.ident))
label <- as.data.frame(label)
label[grep("blasto",label[,1]),2] <- "blasto"
h3k4me1@meta.data$group.ident <- as.factor(label[,2])
DimPlot(h3k4me1, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off()

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me1/h3k4me1_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(h3k4me1)]
h3k4me1[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k4me1) <- "ACTIVITY"
h3k4me1 <- FindVariableFeatures(h3k4me1, nfeatures = 10000)
h3k4me1 <- NormalizeData(h3k4me1)
h3k4me1 <- ScaleData(h3k4me1)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = blasto_pseudo,
                                        query = h3k4me1,
                                        features = VariableFeatures(object = blasto_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 10, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)

#transfer ID
blasto_pseudo@meta.data$ID <- colnames(blasto_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = blasto_pseudo$ID,
                                   weight.reduction = h3k4me1[["lsi"]], dims = 2:ncol(h3k4me1[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k4me1_pre1 <- ID.predict
length(unique(h3k4me1_pre1$predicted.id))
#h3k4me1_pre1 <- h3k4me1_pre1[which(h3k4me1_pre1$prediction.score.max >0.2),]
write.table(h3k4me1_pre1, "../h3k4me1/h3k4me1_blasto_pseudo.txt", quote = F, sep = " ")

saveRDS(h3k4me1, "../h3k4me1/h3k4me1_blasto_integration.rds")
saveRDS(transfer.anchors, "../h3k4me1/h3k4me1_blasto_transfer_anchors.rds")

######################################### H3K27ac #######################################
setwd("../h3k27ac/")
set.seed(2)
h3k27ac_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/h3k27ac.5kb.best.countMatrix.txt", sep =" ") 
h3k27ac_counts <- h3k27ac_counts[,grep("blasto|E3.5|E3.0",colnames(h3k27ac_counts))] 
h3k27ac <- CreateSeuratObject(counts = h3k27ac_counts,assay = 'peaks',project = 'h3k27ac',min.cells = 1)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = 'q0')
h3k27ac <- RunSVD(object = h3k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k27ac/h3k27ac_allStage_peaks_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,grep("blasto|E3.5|E3.0", colnames(geneActivity))]

dat <- data.frame(col=colnames(geneActivity), min="min")
dat_1 <- dat[grep("talib2238", dat$col),]
dat_2 <- dat[grep("talib2238", dat$col, invert = T),]
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
transfer.anchors <- FindTransferAnchors(reference = blasto_pseudo,
                                        query = h3k27ac,
                                        features = VariableFeatures(object = blasto_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 10, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)

#transfer ID
blasto_pseudo@meta.data$ID <- colnames(blasto_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = blasto_pseudo$ID,
                                   weight.reduction = h3k27ac[["lsi"]], dims = 2:ncol(h3k27ac[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k27ac_pre1 <- ID.predict
length(unique(h3k27ac_pre1$predicted.id))
#h3k27ac_pre1 <- h3k27ac_pre1[which(h3k27ac_pre1$prediction.score.max >0.3),]
write.table(h3k27ac_pre1, "./h3k27ac_blasto_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k27ac, "h3k27ac_blasto_integration.rds")
saveRDS(transfer.anchors, "h3k27ac_blasto_transfer_anchors.rds")


############################################### H3K36me3 #########################################
#rm(list = ls())
setwd("../h3k36me3/")
set.seed(3)
h3k36me3 <- readRDS("./h3k36me3_5kb.rds")
h3k36me3 <- subset(h3k36me3,cells = colnames(h3k36me3[,grep("E3.0|E3.5", h3k36me3@meta.data$group.ident)]))
h3k36me3 <- RunTFIDF(h3k36me3)
h3k36me3 <- FindTopFeatures(h3k36me3, min.cutoff = 'q0')
h3k36me3 <- RunSVD(object = h3k36me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
h3k36me3 <- RunUMAP(object = h3k36me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
h3k36me3 <- FindNeighbors(object = h3k36me3,reduction = 'lsi',dims = c(2:15))
h3k36me3 <- FindClusters(object = h3k36me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
pdf("h3k36me3_blasto_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k36me3/h3k36me3_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,grep("E3.0|E3.5", colnames(geneActivity))]
geneActivity <- geneActivity[,colnames(h3k36me3)]
h3k36me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k36me3) <- "ACTIVITY"
h3k36me3 <- FindVariableFeatures(h3k36me3, nfeatures = 10000)
h3k36me3 <- NormalizeData(h3k36me3)
h3k36me3 <- ScaleData(h3k36me3)

#RNA <- readRDS("../RNA_blasto.rds")
#####integration#####
transfer.anchors <- FindTransferAnchors(reference = blasto_pseudo,
                                        query = h3k36me3,
                                        features = VariableFeatures(object = blasto_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 10, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
blasto_pseudo@meta.data$ID <- colnames(blasto_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = blasto_pseudo$ID,
                                   weight.reduction = h3k36me3[["lsi"]], dims = 2:ncol(h3k36me3[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k36me3_pre1 <- ID.predict
length(unique(h3k36me3_pre1$predicted.id))
#h3k36me3_pre1 <- h3k36me3_pre1[which(h3k36me3_pre1$prediction.score.max >0.2),]
write.table(h3k36me3_pre1, "./h3k36me3_blasto_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k36me3, "h3k36me3_blasto_integration.rds")
saveRDS(transfer.anchors, "h3k36me3_blasto_transfer_anchors.rds")

########################################### H3K27me3 #############################################
##借助coTACIT数据
setwd("../cotacit_data/")
cotacit <- readRDS("./cotacit_h3k27ac.rds")
cotacit_blasto <- subset(cotacit,cells = colnames(cotacit[,grep("blasto", cotacit@meta.data$group.ident)]))
geneActivity <- read.csv("./cotacit_h3k27ac_allstage_tss5kb_ActivityMatrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(cotacit_blasto)]
cotacit_blasto[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(cotacit_blasto) <- "ACTIVITY"
cotacit_blasto <- FindVariableFeatures(cotacit_blasto, nfeatures = 10000)
cotacit_blasto <- NormalizeData(cotacit_blasto)
cotacit_blasto <- ScaleData(cotacit_blasto)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = blasto_pseudo,
                                        query = cotacit_blasto,
                                        features = VariableFeatures(object = blasto_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 10, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
blasto_pseudo@meta.data$ID <- colnames(blasto_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = blasto_pseudo$ID,
                                   weight.reduction = cotacit_blasto[["lsi"]], dims = 2:ncol(cotacit_blasto[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_pre1 <- ID.predict
length(unique(cotacit_pre1$predicted.id))
#cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]
write.table(cotacit_pre1, "./cotacit_blasto_h3k27ac_pseudo.txt", quote = F, sep = " ")
saveRDS(cotacit_blasto, "cotacit_blasto_h3k27ac_integration.rds")
saveRDS(transfer.anchors, "cotacit_blasto_h3k27ac_transfer_anchors.rds")


####################################### all #######################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_blasto_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_blasto_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_blasto_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_blasto_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_blasto_h3k27ac_pseudo.txt", sep = " ")

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

write.table(h3k4me1_pre1, "./h3k4me1_blasto_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k4me3_pre1, "./h3k4me3_blasto_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k27ac_pre1, "./h3k27ac_blasto_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k36me3_pre1, "./h3k36me3_blasto_pseudo_select.txt", quote = F, sep = " ")
write.table(cotacit_pre1, "./cotacit_blasto_pseudo_select.txt", quote = F, sep = " ")
##由于想尽可能多地保留整合后的细胞，所以整合到同一个RNA细胞的多个细胞也分开并保留，因此在predicted.id后面还增加了一个编号，这一步在excel中操作(按照predicted.id和score值进行sort后编号)。
h3k4me1_pre1 <- read.csv("./h3k4me1_blasto_pseudo_select.txt", sep = "\t")
h3k4me3_pre1 <- read.csv("./h3k4me3_blasto_pseudo_select.txt", sep = "\t")
h3k27ac_pre1 <- read.csv("./h3k27ac_blasto_pseudo_select.txt", sep = "\t")
h3k36me3_pre1 <- read.csv("./h3k36me3_blasto_pseudo_select.txt", sep = "\t")
cotacit_pre1 <- read.csv("./cotacit_blasto_pseudo_select.txt", sep = "\t")

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))
inte
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1$predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1$predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1$predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1$predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$predicted.id %in% inte),]

h3k4me1.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/", h3k4me1_pre1$X, " ./blasto/", h3k4me1_pre1$predicted.id, "_inte.bam")
write.table(h3k4me1.sh, "./h3k4me1/h3k4me1_blasto_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", h3k4me3_pre1$X, " ./blasto/", h3k4me3_pre1$predicted.id, "_inte.bam")
write.table(h3k4me3.sh, "./h3k4me3/h3k4me3_blasto_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k27ac_pre1_1 <- h3k27ac_pre1[grep("talib2238", h3k27ac_pre1$X),]
h3k27ac_pre1_2 <- h3k27ac_pre1[grep("talib2238", h3k27ac_pre1$X, invert = T),]
h3k27ac_pre1_1 <- separate(h3k27ac_pre1_1, X, c("x1", "x2", "x3", "x4", "x5"), sep = "\\.")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, id, c("x1", "x2", "x3"), sep = "-")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, X, c("id", "x4", "x5"), sep = ".")
h3k27ac_pre1 <- rbind(h3k27ac_pre1_1, h3k27ac_pre1_2)
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", h3k27ac_pre1$X, " ./blasto/", h3k27ac_pre1$predicted.id, "_inte.bam")
write.table(h3k27ac.sh, "./h3k27ac/h3k27ac_blasto_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", h3k36me3_pre1$X, " ./blasto/", h3k36me3_pre1$predicted.id, "_inte.bam")
write.table(h3k36me3.sh, "./h3k36me3/h3k36me3_blasto_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)


############################################# cotacit + tacit h3k37me3 #########################################

############ cotacit h3k27me3 ################
library(Signac)
h3k27me3 <- readRDS("./h3k27me3_all stage_peaks.rds")
cotacit_h3k27me3 <- readRDS("./cotacit_data/cotacit_h3k27me3_peaks.rds")

h3k27me3$dataset <- "tacit"
cotacit_h3k27me3$dataset <- "cotacit"

h3k27me3_blasto <- subset(h3k27me3,cells = colnames(h3k27me3)[grep("E3.5|E3.0|blasto", colnames(h3k27me3))])
cotacit_h3k27me3_blasto <- subset(cotacit_h3k27me3,cells = colnames(cotacit_h3k27me3)[grep("blasto", cotacit@meta.data$group.ident)])
merge <- merge(h3k27me3_blasto, cotacit_h3k27me3_blasto)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_blasto.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k27me3_blasto, cotacit_h3k27me3_blasto),
  anchor.features = rownames(h3k27me3_blasto),
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

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_integration_blasto.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k27me3_blasto,
  query = cotacit_h3k27me3_blasto,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k27me3_blasto@meta.data$ID <- colnames(h3k27me3_blasto)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k27me3_blasto$ID,
                                   weight.reduction = cotacit_h3k27me3_blasto[["lsi"]], dims = 2:ncol(cotacit_h3k27me3_blasto[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k27me3_pre1 <- ID.predict
length(unique(cotacit_h3k27me3_pre1$predicted.id))

write.table(cotacit_h3k27me3_pre1, "./cotacit_data/h3k27me3_blasto_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k27me3_blasto_cotacit_tacit_inte_transfer_anchors.rds")

############################################# cotacit + tacit h3k9me3 #########################################

############ cotacit h3k9me3 ################
cotacit_h3k9me3 <- readRDS("./cotacit_data/cotacit_h3k9me3_10kbbin.rds")
h3k9me3 <- readRDS("./cotacit_data/h3k9me3_10kbbin.rds")
h3k9me3$dataset <- "tacit"
cotacit_h3k9me3$dataset <- "cotacit"

h3k9me3_blasto <- subset(h3k9me3,cells = colnames(h3k9me3)[grep("blasto|E3.5|E3.0", h3k9me3@meta.data$group.ident)])
cotacit_h3k9me3_blasto <- subset(cotacit_h3k9me3,cells = colnames(cotacit_h3k9me3)[grep("blasto|E3.5|E3.0", cotacit_h3k9me3@meta.data$group.ident)])
merge <- merge(h3k9me3_blasto, cotacit_h3k9me3_blasto)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_blasto.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k9me3_blasto, cotacit_h3k9me3_blasto),
  anchor.features = rownames(h3k9me3_blasto),
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

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_integration_blasto.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k9me3_blasto,
  query = cotacit_h3k9me3_blasto,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k9me3_blasto@meta.data$ID <- colnames(h3k9me3_blasto)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k9me3_blasto$ID,
                                   weight.reduction = cotacit_h3k9me3_blasto[["lsi"]], dims = 2:ncol(cotacit_h3k9me3_blasto[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k9me3_pre1 <- ID.predict
length(unique(cotacit_h3k9me3_pre1$predicted.id))

write.table(cotacit_h3k9me3_pre1, "./cotacit_data/h3k9me3_blasto_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k9me3_blasto_cotacit_tacit_inte_transfer_anchors.rds")


cotacit_pre1$ID <- cotacit_pre1$X
cotacit_pre1_1 <- cotacit_pre1[grep("talib23024", cotacit_pre1$X),]
cotacit_pre1_2 <- cotacit_pre1[grep("talib23024", cotacit_pre1$X, invert = T),]

cotacit_pre1_1 <- separate(cotacit_pre1_1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1_1 <- unite(cotacit_pre1_1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1_1$ID2 <- paste0(cotacit_pre1_1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1_1$ID3 <- paste0(cotacit_pre1_1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam") 
cotacit_pre1_1$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1_1$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1_1$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1_1$ID3, rownames(cotacit_h3k9me3_pre1))]
k27me3.sh <- paste0("samtools merge ./h3k27me3_tacit+cotacit/blasto/", cotacit_pre1_1$predicted.id, "_inte.bam ", 
                    "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1_1$ID2,
                    " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1_1$k27_tacit)
k9me3.sh <- paste0("samtools merge ./h3k9me3_tacit+cotacit/blasto/", cotacit_pre1_1$predicted.id, "_inte.bam ", 
                   "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1_1$ID3,
                   " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1_1$k9_tacit)
write.table(k27me3.sh, "./cotacit_data/integration/h3k27me3_blasto_merge_1.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./cotacit_data/integration/h3k9me3_blasto_merge_1.sh", sep = "\t", quote = F, col.names = F, row.names = F)


cotacit_pre1_2$ID2 <- cotacit_pre1_2$ID
cotacit_pre1_2$ID3 <- cotacit_pre1_2$ID
cotacit_pre1_2$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1_2$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1_2$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1_2$ID3, rownames(cotacit_h3k9me3_pre1))]
k27me3.sh <- paste0("samtools merge ./h3k27me3_tacit+cotacit/blasto/", cotacit_pre1_2$predicted.id, "_inte.bam ", 
                    "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1_2$ID2,
                    " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1_2$k27_tacit)
k9me3.sh <- paste0("samtools merge ./h3k9me3_tacit+cotacit/blasto/", cotacit_pre1_2$predicted.id, "_inte.bam ", 
                   "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1_2$ID3,
                   " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1_2$k9_tacit)
write.table(k27me3.sh, "./cotacit_data/integration/h3k27me3_blasto_merge_2.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./cotacit_data/integration/h3k9me3_blasto_merge_2.sh", sep = "\t", quote = F, col.names = F, row.names = F)


################################## 为了用尽可能多的细胞构建chromHMM，基本保留所有整合的TACIT和cotacit细胞 ###########################################
###### 这里和pseudo_4cell_meta.R做的事情类似 #####
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_blasto_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_blasto_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_blasto_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_blasto_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_blasto_h3k27ac_pseudo.txt", sep = " ")

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
write.table(h3k4me1.sh, "./09_chromHMM_blasto/h3k4me1_blasto_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k4me3 ####
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", rownames(h3k4me3_pre1), 
                     " ./", h3k4me3_pre1$predicted.id, "/")
write.table(h3k4me3.sh, "./09_chromHMM_blasto/h3k4me3_blasto_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k27ac ####
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k27ac_pre1$X <- rownames(h3k27ac_pre1)
h3k27ac_pre1_1 <- h3k27ac_pre1[grep("2238", rownames(h3k27ac_pre1)),]
h3k27ac_pre1_2 <- h3k27ac_pre1[grep("2238", rownames(h3k27ac_pre1), invert = T),]

h3k27ac_pre1_1 <- separate(h3k27ac_pre1_1, X, c("x1", "x2", "x3", "x4", "x5"), sep = "\\.")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, id, c("x1", "x2", "x3"), sep = "-")
h3k27ac_pre1_1 <- unite(h3k27ac_pre1_1, X, c("id", "x4", "x5"), sep = ".")
h3k27ac_pre1 <- rbind(h3k27ac_pre1_1, h3k27ac_pre1_2)
rownames(h3k27ac_pre1) <- h3k27ac_pre1$X
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", rownames(h3k27ac_pre1), 
                     " ./", h3k27ac_pre1$predicted.id, "/")
write.table(h3k27ac.sh, "./09_chromHMM_blasto/h3k27ac_blasto_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k36me3 ####
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", rownames(h3k36me3_pre1), 
                      " ./", h3k36me3_pre1$predicted.id, "/")
write.table(h3k36me3.sh, "./09_chromHMM_blasto/h3k36me3_blasto_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)


#### cotacit ####
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_blasto_h3k27ac_pseudo.txt", sep = " ")
cotacit_h3k27me3_pre1 <- read.csv("./cotacit_data/h3k27me3_blasto_cotacit_tacit_inte.txt", sep = " ")
cotacit_h3k9me3_pre1 <- read.csv("./cotacit_data/h3k9me3_blasto_cotacit_tacit_inte.txt", sep = " ")

cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$predicted.id %in% inte),]
cotacit_pre1$ID <- rownames(cotacit_pre1)
cotacit_pre1_1 <- cotacit_pre1[grep("talib23024", rownames(cotacit_pre1)),]
cotacit_pre1_2 <- cotacit_pre1[grep("talib23024", rownames(cotacit_pre1), invert = T),]

cotacit_pre1_1 <- separate(cotacit_pre1_1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1_1 <- unite(cotacit_pre1_1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1_1$ID2 <- paste0(cotacit_pre1_1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1_1$ID3 <- paste0(cotacit_pre1_1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam") 
cotacit_pre1_1$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1_1$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1_1$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1_1$ID3, rownames(cotacit_h3k9me3_pre1))]

k27me3_tacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1_1$k27_tacit, 
                          " ./", cotacit_pre1_1$predicted.id, "/")
k27me3_cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1_1$ID2, 
                            " ./", cotacit_pre1_1$predicted.id, "/")
k9me3_tacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1_1$k9_tacit, 
                         " ./", cotacit_pre1_1$predicted.id, "/")
k9me3_cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1_1$ID3, 
                           " ./", cotacit_pre1_1$predicted.id, "/")
k27me3.sh <- c(k27me3_tacit.sh, k27me3_cotacit.sh)
k9me3.sh <- c(k9me3_tacit.sh, k9me3_cotacit.sh)
write.table(k27me3.sh, "./09_chromHMM_blasto/h3k27me3_blasto_meta_cp_1.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./09_chromHMM_blasto/h3k9me3_blasto_meta_cp_1.sh", sep = "\t", quote = F, col.names = F, row.names = F)

cotacit_pre1_2$ID2 <- cotacit_pre1_2$ID
cotacit_pre1_2$ID3 <- cotacit_pre1_2$ID
cotacit_pre1_2$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1_2$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1_2$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1_2$ID3, rownames(cotacit_h3k9me3_pre1))]
k27me3_tacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1_2$k27_tacit, 
                          " ./", cotacit_pre1_2$predicted.id, "/")
k27me3_cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/coTarget_blasto/h3k27me3/allbam/", cotacit_pre1_2$ID2, 
                            " ./", cotacit_pre1_2$predicted.id, "/")
k9me3_tacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1_2$k9_tacit, 
                         " ./", cotacit_pre1_2$predicted.id, "/")
k9me3_cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/coTarget_blasto/h3k9me3/allbam/", cotacit_pre1_2$ID3, 
                           " ./", cotacit_pre1_2$predicted.id, "/")
k27me3.sh <- c(k27me3_tacit.sh, k27me3_cotacit.sh)
k9me3.sh <- c(k9me3_tacit.sh, k9me3_cotacit.sh)
write.table(k27me3.sh, "./09_chromHMM_blasto/h3k27me3_blasto_meta_cp_2.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./09_chromHMM_blasto/h3k9me3_blasto_meta_cp_2.sh", sep = "\t", quote = F, col.names = F, row.names = F)



