########################scRNA-seq datasets comparision #################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA")

library(Seurat)#version 4.3.0
library(SeuratObject)#version 4.1.3
library(dplyr)
set.seed(123)

RNA <- CreateSeuratObject(counts = RNA_counts, project = "RNA", min.cells = 30, min.features = 3000)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = "^MT-")
pdf("all_RNA_no_outlier.pdf")
VlnPlot(RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
plot1 <- FeatureScatter(RNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
plot2
RNA <- NormalizeData(RNA)
RNA <- FindVariableFeatures(RNA, nfeatures = 5000)
RNA <- ScaleData(RNA, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
RNA <- RunPCA(RNA, npcs = 20)
RNA <- RunUMAP(RNA, dims = 1:20, reduction = 'pca',n.neighbors = 100L)
RNA <- FindNeighbors(RNA, dims = 1:20, reduction = "pca")
RNA <- FindClusters(RNA, resolution = 0.5)
DimPlot(RNA, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(RNA, reduction = "pca", label = TRUE, pt.size = 2)
dev.off()

label <- rownames(as.data.frame(RNA@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "late2cell"
label[grep("early2cell",label[,1]),2] <- "early2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("morula",label[,1]),2] <- "Morula"
label[grep("earlyblast",label[,1]),2] <- "earlyblast"
label[grep("lateblast",label[,1]),2] <- "lateblast"

label[,3] <- RNA@meta.data$seurat_clusters
RNA@meta.data$group.ident <- as.factor(label[,2])
RNA@meta.data$ID <- as.factor(label[,1])

pdf("RNA_stage.pdf")
DimPlot(RNA, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(RNA, reduction = "pca", label = TRUE, pt.size = 2)
DimPlot(object = RNA, reduction = "pca",group.by = "group.ident" ,label = TRUE, pt.size = 2)
dev.off()

# find markers
RNA.markers <- FindAllMarkers(RNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
RNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50
pdf("stage_marker_heatmap.pdf")
DoHeatmap(RNA, features = top50$gene) + NoLegend()
dev.off()
saveRDS(RNA, file = "./all_RNA.rds")

#### reference scRNA-seq ####
ref <- read.csv("./scRNA-seq_2014/GSE45719_scRNA.csv")
ref <- ref[!duplicated(ref$Gene_symbol),]
rownames(ref) <- ref[,1]
ref <- ref[,-1]
ref <- as.matrix(ref)

fpkmToTpm <- function(fpkm) 
{ 
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6)) 
}
ref_tpm <- apply(ref,2,fpkmToTpm)

ref <- CreateSeuratObject(counts = ref_tpm, project = "ref", min.cells = 2)
ref[["percent.mt"]] <- PercentageFeatureSet(ref, pattern = "^MT-")
pdf("all_reference.pdf")
VlnPlot(ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
plot1 <- FeatureScatter(ref, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ref, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
plot2
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref, nfeatures = 5000)
ref <- ScaleData(ref, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
ref <- RunPCA(ref, npcs = 20)
ref <- RunUMAP(ref, dims = 1:20, reduction = 'pca',n.neighbors = 5L)
ref <- FindNeighbors(ref, dims = 1:20, reduction = "pca")
ref <- FindClusters(ref, resolution = 0.5)
DimPlot(ref, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(ref, reduction = "pca", label = TRUE, pt.size = 2)
dev.off()

pdf("ref_stage.pdf")
DimPlot(ref, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = ref, reduction = "umap",group.by = "orig.ident" ,label = TRUE, pt.size = 2)
DimPlot(ref, reduction = "pca", label = TRUE, pt.size = 2)
DimPlot(object = ref, reduction = "pca",group.by = "orig.ident" ,label = TRUE, pt.size = 2)
dev.off()

saveRDS(ref, file = "./reference.rds")

###################################### integration ##############################################
merge <- merge(RNA, y = ref, add.cell.ids = c("smart-seq3", "smart-seq2"), project = "embryo")
embryo.list <- SplitObject(merge, split.by = "dataset")
embryo.list

# normalization
embryo.list <- lapply(X = embryo.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})
# integration
features <- SelectIntegrationFeatures(object.list = embryo.list)
embryo.anchors <- FindIntegrationAnchors(object.list = embryo.list, anchor.features = features)
embryo.combined <- IntegrateData(anchorset = embryo.anchors)

DefaultAssay(embryo.combined) <- "integrated"
embryo.combined <- ScaleData(embryo.combined, verbose = FALSE)
embryo.combined <- RunPCA(embryo.combined, npcs = 30, verbose = FALSE)
embryo.combined <- RunUMAP(embryo.combined, reduction = "pca", dims = 1:30)
embryo.combined <- FindNeighbors(embryo.combined, reduction = "pca", dims = 1:30)
embryo.combined <- FindClusters(embryo.combined, resolution = 0.5)

pdf("integrate.pdf")
DimPlot(embryo.combined, reduction = "umap", group.by = "dataset")
DimPlot(embryo.combined, reduction = "umap", group.by = "group.ident")
DimPlot(embryo.combined, reduction = "umap", label = TRUE, repel = TRUE)
dev.off()

saveRDS(embryo.combined, file = "./embryo.combined.rds")

################ SCTransfermation #####################
# library(sctransform)
# 
# embryo.list <- SplitObject(merge, split.by = "dataset")
# embryo.list <- lapply(X = embryo.list, FUN = function(x) {
#   x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent.mt")
#   x <- SCTransform(x, vars.to.regress = "percent.mt", verbose = FALSE)
# })
# 
# features <- SelectIntegrationFeatures(object.list = embryo.list, nfeatures = 3000)
# embryo.list <- PrepSCTIntegration(object.list = embryo.list, anchor.features = features)
# embryo.anchors <- FindIntegrationAnchors(object.list = embryo.list, normalization.method = "SCT",
#                                          anchor.features = features)
# embryo.combined.sct <- IntegrateData(anchorset = embryo.anchors, normalization.method = "SCT")
# embryo.combined.sct <- RunPCA(embryo.combined.sct, verbose = FALSE)
# embryo.combined.sct <- RunUMAP(embryo.combined.sct, reduction = "pca", dims = 1:30)
# 
# pdf("integrate.pdf")
# DimPlot(embryo.combined.sct, reduction = "umap", group.by = "dataset")
# DimPlot(embryo.combined.sct, reduction = "umap", group.by = "group.ident")
# DimPlot(embryo.combined.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
#               repel = TRUE)
# dev.off()



#################################### ZGA gene expression in early-2cells ############################
RNA <- readRDS("./all_RNA_no_outlier.rds")
RNA_count <- as.data.frame(RNA@assays[["RNA"]]@counts)
RNA_count <- RNA_count[which(rowSums(RNA_count) >0),]
toti <- read.csv("../totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("../pluripotency.markers.csv", sep = ",", header = F)
ZGA <- read.xlsx("./MajorZGA_genes.xlsx")
ZGA <- na.omit(ZGA$Gene)
RNA_ZGA <- RNA_count[which(rownames(RNA_count) %in% ZGA),]
RNA_toti <- RNA_count[which(rownames(RNA_count) %in% toti$V1),]
RNA_pluri <- RNA_count[which(rownames(RNA_count) %in% pluri$V1),]
RNA_ZGA <- colSums(RNA_ZGA)/nrow(RNA_ZGA)
RNA_toti <- colSums(RNA_toti)/nrow(RNA_toti)
RNA_pluri <- colSums(RNA_pluri)/nrow(RNA_pluri)

label <- data.frame(label=colnames(RNA_count))
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "late2cell"
label[grep("early2cell",label[,1]),2] <- "early2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("morula",label[,1]),2] <- "Morula"
label[grep("earlyblast",label[,1]),2] <- "earlyblast"
label[grep("lateblast",label[,1]),2] <- "lateblast"
label$V2 <- factor(label$V2, levels=c("zygote", "early2cell", "late2cell", "4cell", "8cell", "Morula", "earlyblast", "lateblast"))

label$ZGA <- RNA_ZGA
label$toti <- RNA_toti
label$pluri <- RNA_pluri
pdf("./RNA_allstage_zga_average.pdf")
ggplot(label,mapping = aes(x = V2, y = ZGA))+
  geom_violin(aes(fill =V2), trim = T) + 
  geom_boxplot(width = 0.05)+
  ggtitle("ZGA genes")+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0","#AC73EF", "#2A2AF4"))+
  theme(legend.position = "none")+
  theme_bw()
dev.off()

################################################# generate RNA pseudo cells ###############################################

############################################# RNA ########################################################
RNA <- readRDS("./all_RNA.rds")

#### monocle for defining pseudo cells #####
library(monocle3)
set.seed(1234)
data <- GetAssayData(RNA, assay = 'RNA', slot = 'counts')
label <- as.data.frame(colnames(RNA))
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "late2cell"
label[grep("early2cell",label[,1]),2] <- "early2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("morula",label[,1]),2] <- "Morula"
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
pdf("RNA_trajectory.pdf")
p1
p2

## trajectory
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "stage", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE,  graph_label_size=3, cell_size = 1.2)
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, label_leaves=TRUE, label_branch_points=TRUE,  graph_label_size=3, cell_size = 1.2)
dev.off()

pseudo_order <- sort(cds@principal_graph_aux@listData$UMAP$pseudotime)
pseudo_order <- as.data.frame(pseudo_order)
saveRDS(cds, "./RNA_trajectory.rds")

RNA_counts <- as.matrix(RNA@assays[["RNA"]]@counts)
RNA_counts <- RNA_counts[,rownames(pseudo_order)]

a <- rep(1:80, each=5)
write.table(a,"./scRNA/tmp.txt", quote = F, sep = "\t")
RNA_counts_zygote <- RNA_counts[,grep("zygote", colnames(RNA_counts))]
RNA_counts_2cell <- RNA_counts[,grep("2cell", colnames(RNA_counts))]
RNA_counts_4cell <- RNA_counts[,grep("4cell", colnames(RNA_counts))]
RNA_counts_8cell <- RNA_counts[,grep("8cell", colnames(RNA_counts))]
RNA_counts_morula <- RNA_counts[,grep("morula", colnames(RNA_counts))]
RNA_counts_blasto <- RNA_counts[,grep("blasto", colnames(RNA_counts))]

# pesudo-cells, merge every 5 cells, for each stage one by one
dim(RNA_counts_zygote)#16619   136
dat3 <- matrix(1:315761,ncol=19)
id <- seq(from=1,to=95,by=5)
n=1
for (i in id) {
  m=i+4
  dat3[,n] <- rowSums(RNA_counts_zygote[,i:m])
  n=n+1
}
dat3 <- as.data.frame(dat3)
rownames(dat3) <- rownames(RNA_counts_zygote)
colnames(dat3) <- paste0("zygote_pseudo_", 1:ncol(dat3))
pseudo_zygote <- dat3 

RNA_pseudo <- cbind(pseudo_zygote, pseudo_2cell, pseudo_4cell, pseudo_8cell, pseudo_morula, pseudo_blasto)

RNA_pseudo <- CreateSeuratObject(counts = RNA_pseudo, project = "RNA", min.cells = 3, min.features = 3000)
RNA_pseudo[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo, pattern = "^MT-")
pdf("RNA_pseudo.pdf")
VlnPlot(RNA_pseudo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
plot1 <- FeatureScatter(RNA_pseudo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(RNA_pseudo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
plot2
RNA_pseudo <- NormalizeData(RNA_pseudo)
RNA_pseudo <- FindVariableFeatures(RNA_pseudo, nfeatures = 8000)
RNA_pseudo <- ScaleData(RNA_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_pseudo <- RunPCA(RNA_pseudo, npcs = 15)
RNA_pseudo <- RunUMAP(RNA_pseudo, dims = 1:15, reduction = 'pca',n.neighbors = 5L)
RNA_pseudo <- FindNeighbors(RNA_pseudo, dims = 1:15, reduction = "pca")
RNA_pseudo <- FindClusters(RNA_pseudo, resolution = 0.5)
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(RNA_pseudo, reduction = "pca", label = TRUE, pt.size = 2)
dev.off()

label <- rownames(as.data.frame(RNA_pseudo@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "late2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("morula",label[,1]),2] <- "Morula"
label[grep("blast",label[,1]),2] <- "blast"

label[,3] <- RNA_pseudo@meta.data$seurat_clusters
RNA_pseudo@meta.data$group.ident <- as.factor(label[,2])
RNA_pseudo@meta.data$ID <- as.factor(label[,1])

pdf("RNA_pseudo_stage.pdf")
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(RNA_pseudo, reduction = "pca", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, reduction = "pca",group.by = "group.ident" ,label = TRUE, pt.size = 2)
dev.off()

saveRDS(RNA_pseudo, "./scRNA/RNA_pseudo.rds")
write.table(RNA_pseudo, "./scRNA/RNA_speudo_count.txt", quote = F, sep = "\t")





