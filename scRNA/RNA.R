########################scRNA-seq datasets comparision #################################

setwd("./scRNA")
############## Seurat ####################
library(Seurat)#version 4.3.0
library(SeuratObject)#version 4.1.3
library(dplyr)
set.seed(123)

RNA <- CreateSeuratObject(counts = RNA_tpm, project = "RNA", min.cells = 30, min.features = 3000)
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

pdf("RNA_stage_no_outlier.pdf")
DimPlot(RNA, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(RNA, reduction = "pca", label = TRUE, pt.size = 2)
DimPlot(object = RNA, reduction = "pca",group.by = "group.ident" ,label = TRUE, pt.size = 2)
dev.off()

RNA.markers <- FindAllMarkers(RNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
RNA.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50
pdf("stage_marker_heatmap_no_outlier.pdf")
DoHeatmap(RNA, features = top50$gene) + NoLegend()
dev.off()
saveRDS(RNA, file = "./all_RNA.rds")

#### reference scRNA-seq ####
ref <- read.csv("./GSE45719_scRNA.csv")
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

################################################################################################
####################################### combind ###############################################
row <- intersect(rownames(ref_tpm), rownames(RNA_no_outlier))
ref_tpm_inte <- ref_tpm[row,]
RNA_no_outlier_inte <- RNA_no_outlier[row,]
all_data <- cbind(RNA_no_outlier_inte, ref_tpm_inte) 

merge <- CreateSeuratObject(counts = all_data, project = "merge", min.cells = 20)
merge[["percent.mt"]] <- PercentageFeatureSet(merge, pattern = "^MT-")
pdf("merge_RNA.pdf")
VlnPlot(merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)
plot1 <- FeatureScatter(merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
plot2
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge, nfeatures = 5000)
merge <- ScaleData(merge, vars.to.regress = c("nFeature_RNA", "percent.mito","S.Score", "G2M.Score"))
merge <- RunPCA(merge, npcs = 30)
merge <- RunUMAP(merge, dims = 1:30, reduction = 'pca',n.neighbors = 100L)
merge <- FindNeighbors(merge, dims = 1:30, reduction = "pca")
merge <- FindClusters(merge, resolution = 0.5)
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(merge, reduction = "pca", label = TRUE, pt.size = 2)
dev.off()

label <- rownames(as.data.frame(merge@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "late2cell"
label[grep("early2cell",label[,1]),2] <- "early2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("morula|16cell",label[,1]),2] <- "Morula"
label[grep("earlyblast",label[,1]),2] <- "earlyblast"
label[grep("midblast",label[,1]),2] <- "midblast"
label[grep("lateblast",label[,1]),2] <- "lateblast"

label[,3] <- merge@meta.data$seurat_clusters
label[grep("CP",label[,1]),4] <- "this study"
label[grep("CP",label[,1], invert = T),4] <- "reference"

merge@meta.data$group.ident <- as.factor(label[,2])
merge@meta.data$ID <- as.factor(label[,1])
merge@meta.data$dataset <- as.factor(label[,4])

pdf("merge_stage.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "orig.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
DimPlot(merge, reduction = "pca", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "pca",group.by = "orig.ident" ,label = TRUE, pt.size = 2)
dev.off()

####################################### integration ##############################################
embryo.list <- SplitObject(merge, split.by = "dataset")
embryo.list

embryo.list <- lapply(X = embryo.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
})

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

library(openxlsx)
#################################### ZGA gene expression in early-2cells ############################
RNA <- readRDS("./all_RNA.rds")
RNA_count <- as.data.frame(RNA@assays[["RNA"]]@counts)
RNA_count <- RNA_count[which(rowSums(RNA_count) >0),]
toti <- read.csv("./totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("./pluripotency.markers.csv", sep = ",", header = F)
ZGA <- read.xlsx("./MajorZGA_genes.xlsx")
ZGA <- na.omit(ZGA$Gene)
RNA_ZGA <- RNA_count[which(rownames(RNA_count) %in% ZGA),]
RNA_toti <- RNA_count[which(rownames(RNA_count) %in% toti$V1),]
RNA_pluri <- RNA_count[which(rownames(RNA_count) %in% pluri$V1),]
RNA_ZGA <- colSums(RNA_ZGA)/620
RNA_toti <- colSums(RNA_toti)/215
RNA_pluri <- colSums(RNA_pluri)/38

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








