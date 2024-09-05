

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

######################## RNA #########################
RNA <- readRDS("./scRNA/all_RNA.rds")
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

################################ generate ICM and TE pseudo cells#########################
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
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('cds.umap')
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(RNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('int.umap')
pdf("./scRNA/RNA_ICM_trajectory.pdf")
p1
p2

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

cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('cds.umap')

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(RNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('int.umap')
pdf("./scRNA/RNA_TE_trajectory.pdf")
p1
p2

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

############################################# H3K4me3 (as an example) ####################################################
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

####################################### all (only retain cells with prediction.score.max >0.2 across all modalities) #######################################################
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
