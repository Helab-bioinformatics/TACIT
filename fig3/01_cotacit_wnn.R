setwd("./early2cell/WNN")
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

############# H3K4me3 #############
pathToBams1 <- c('../h3k4me3/rename/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
regions <- "./mm10.5K.windows.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3')
h3k4me3_counts <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
h3k4me3 <- CreateSeuratObject(
  counts = h3k4me3_counts,
  assay = 'peaks',
  project = 'h3k4me3',
  min.cells = 1
)
h3k4me3 <- RunTFIDF(h3k4me3)
h3k4me3 <- FindTopFeatures(h3k4me3, min.cutoff = 'q0')
h3k4me3 <- RunSVD(object = h3k4me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("h3k4me3_early2cell.pdf")
DepthCor(h3k4me3)
h3k4me3 <- RunUMAP(
  object = h3k4me3,
  reduction = 'lsi',
  dims = c(2:20)
)
h3k4me3 <- FindNeighbors(
  object = h3k4me3,
  reduction = 'lsi',
  dims = c(2:20)
)
h3k4me3 <- FindClusters(
  object = h3k4me3,
  algorithm = 3,
  resolution = 0.8,
  verbose = FALSE
)
DimPlot(object = h3k4me3, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

############### H3K27ac #################
pathToBams2 <- c('../h3k27ac/rename/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
regions <- "./mm10.2K.windows.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles2, regions, paired = T, project.name='early2cell_h3k27ac')
h3k27ac_counts <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
h3k27ac <- CreateSeuratObject(
  counts = h3k27ac_counts,
  assay = 'peaks',
  project = 'h3k27ac',
  min.cells = 1
)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = 'q0')
h3k27ac <- RunSVD(object = h3k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("h3k27ac_early2cell.pdf")
DepthCor(h3k27ac)
h3k27ac <- RunUMAP(object = h3k27ac, reduction = 'lsi', dims = c(2:20))
h3k27ac <- FindNeighbors(object = h3k27ac, reduction = 'lsi',dims = c(2:20))
h3k27ac <- FindClusters(object = h3k27ac, algorithm = 3, resolution = 0.8, verbose = FALSE)
DimPlot(object = h3k27ac, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

##################### H3K27me3 #####################
pathToBams3 <- c('../h3k27me3/rename/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
regions <- "./mm10.5K.windows.bed"
cisTopicObject_h3k27me3 <- createcisTopicObjectFromBAM(bamFiles3, regions, paired = T, project.name='early2cell_h3k27me3')
h3k27me3_counts <- as.matrix(cisTopicObject_h3k27me3@count.matrix)
h3k27me3 <- CreateSeuratObject(
  counts = h3k27me3_counts,
  assay = 'peaks',
  project = 'h3k27me3',
  min.cells = 1
)
h3k27me3 <- RunTFIDF(h3k27me3)
h3k27me3 <- FindTopFeatures(h3k27me3, min.cutoff = 'q0')
h3k27me3 <- RunSVD(object = h3k27me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("h3k27me3_early2cell.pdf")
DepthCor(h3k27me3)
h3k27me3 <- RunUMAP(object = h3k27me3, reduction = 'lsi', dims = c(2:15))
h3k27me3 <- FindNeighbors(object = h3k27me3,reduction = 'lsi',dims = c(2:15))
h3k27me3 <- FindClusters(object = h3k27me3, algorithm = 3, resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k27me3, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

################################ WNN Merge #######################################
coTarget <- CreateSeuratObject(
  counts = h3k4me3_counts,
  assay = 'h3k4me3',
  project = 'h3k4me3',
  min.cells = 1
)
DefaultAssay(coTarget) <- "h3k4me3"
coTarget <- RunTFIDF(coTarget)
coTarget <- FindTopFeatures(coTarget, min.cutoff = 'q0')
coTarget <- RunSVD(object = coTarget, assay = 'h3k4me3', reduction.key = 'h3k4me3lsi_', reduction.name = 'h3k4me3_lsi')
coTarget <- RunUMAP(object = coTarget, reduction = 'h3k4me3_lsi', dims = c(2:20),reduction.name = "umap.h3k4me3", reduction.key = "h3k4me3UMAP_")
coTarget <- FindNeighbors(object = coTarget, reduction = 'h3k4me3_lsi',dims = c(2:20))
coTarget <- FindClusters(object = coTarget,algorithm = 3, resolution = 0.8, verbose = FALSE)

coTarget[["h3k27ac"]] <- h3k27ac@assays$peaks
coTarget[["h3k27me3"]] <- h3k27me3@assays$peaks

DefaultAssay(coTarget) <- "h3k27ac"
coTarget <- RunTFIDF(coTarget)
coTarget <- FindTopFeatures(coTarget, min.cutoff = 'q0')
coTarget <- RunSVD(object = coTarget, assay = 'h3k27ac', reduction.key = 'h3k27aclsi_', reduction.name = 'h3k27ac_lsi')
coTarget <- RunUMAP(object = coTarget, reduction = 'h3k27ac_lsi', dims = c(2:20), reduction.name = "umap.h3k27ac", reduction.key = "h3k27acUMAP_")
coTarget <- FindNeighbors(object = coTarget, reduction = 'h3k27ac_lsi',dims = c(2:20))
coTarget <- FindClusters(object = coTarget,algorithm = 3, resolution = 0.8, verbose = FALSE)

DefaultAssay(coTarget) <- "h3k27me3"
coTarget <- RunTFIDF(coTarget)
coTarget <- FindTopFeatures(coTarget, min.cutoff = 'q0')
coTarget <- RunSVD(object = coTarget, assay = 'h3k27me3', reduction.key = 'h3k27me3lsi_', reduction.name = 'h3k27me3_lsi')
coTarget <- RunUMAP(object = coTarget, reduction = 'h3k27me3_lsi', dims = c(2:15), reduction.name = "umap.h3k27me3", reduction.key = "h3k27me3UMAP_")
coTarget <- FindNeighbors(object = coTarget, reduction = 'h3k27me3_lsi',dims = c(2:15))
coTarget <- FindClusters(object = coTarget,algorithm = 3, resolution = 0.8, verbose = FALSE)

#WNN
coTarget <- FindMultiModalNeighbors(coTarget, reduction.list = list('h3k4me3_lsi', 'h3k27ac_lsi','h3k27me3_lsi'), dims.list = list(2:20, 2:20, 2:20), knn.range = 100)
coTarget <- RunUMAP(coTarget, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
coTarget <- FindClusters(coTarget, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)

coTarget$celltype <- coTarget@meta.data$wsnn_res.0.8
pdf("early2cell_WNN.pdf")
p1 <- DimPlot(coTarget, reduction = "umap.h3k4me3", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k4me3")
p2 <- DimPlot(coTarget, reduction = "umap.h3k27ac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k27ac")
p3 <- DimPlot(coTarget, reduction = "umap.h3k27me3", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k27me3")
p4 <- DimPlot(coTarget, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2+p3+p4+ plot_layout(ncol = 2) & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()
