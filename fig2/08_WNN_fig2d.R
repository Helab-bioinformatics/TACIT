######################################### WNN analysis for interploated single cells (fig. 2d) ########################################
library(cisTopic)
library(ggplot2)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(tidyr)
library(stringr)

############################## h3k4me1 ####################################
pathToBams1 <- c('../2cell/h3k4me1/')
bamFiles1 <- paste(pathToBams1, inte, sep='')
pathToBams2 <- c('../4cell/h3k4me1/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
pathToBams3 <- c('../8cell/h3k4me1/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
pathToBams4 <- c('../../05_morula_inte/h3k4me1/')
bamFiles4 <- paste(pathToBams4, list.files(pathToBams4), sep='')
bamFiles4 <- bamFiles4[grep(".bam", bamFiles4)]
pathToBams5 <- c('../../06_blasto_inte/03_scChromHMM/02_chromHMM_sc/bam/h3k4me1/')
bamFiles5 <- paste(pathToBams5, list.files(pathToBams5), sep='')
pathToBams6 <- c('../zygote/h3k4me1/')
bamFiles6 <- paste(pathToBams6, list.files(pathToBams6), sep='')
bamFiles <- c(bamFiles1, bamFiles2, bamFiles3, bamFiles4, bamFiles5,bamFiles6)

regions <- "./mm10.10K.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, paired = T, project.name='h3k4me1_linked', min.regions = 0)
c <- as.matrix(cisTopicObject@count.matrix) 

h3k4me1 <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
pdf("h3k4me1_10kbbin.pdf")
DepthCor(h3k4me1)
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:15))
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:15))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k4me1@active.ident))
label <- as.data.frame(label)
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("16cell",label[,1]),2] <- "morula"
label[grep("^T5",label[,1]),2] <- "blastocyst"
label[,3] <- h3k4me1@meta.data$seurat_clusters
h3k4me1@meta.data$group.ident <- as.factor(label[,2])
DimPlot(object = h3k4me1, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()


############################## h3k4me3 ####################################
pathToBams1 <- c('../2cell/h3k4me3/')
bamFiles1 <- paste(pathToBams1, inte, sep='')
pathToBams2 <- c('../4cell/h3k4me3/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
pathToBams3 <- c('../8cell/h3k4me3/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
pathToBams4 <- c('../../05_morula_inte/h3k4me3/')
bamFiles4 <- paste(pathToBams4, list.files(pathToBams4), sep='')
bamFiles4 <- bamFiles4[grep(".bam", bamFiles4)]
pathToBams5 <- c('../../06_blasto_inte/03_scChromHMM/02_chromHMM_sc/bam/h3k4me3/')
bamFiles5 <- paste(pathToBams5, list.files(pathToBams5), sep='')
pathToBams6 <- c('../zygote/h3k4me3/')
bamFiles6 <- paste(pathToBams6, list.files(pathToBams6), sep='')
bamFiles <- c(bamFiles1, bamFiles2, bamFiles3, bamFiles4, bamFiles5,bamFiles6)

regions <- "./all_h3k4me3-0.05.5kbMerge.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, paired = T, project.name='h3k4me3_linked', min.regions = 0)
c <- as.matrix(cisTopicObject@count.matrix) 

h3k4me3 <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k4me3',min.cells = 1)
h3k4me3 <- RunTFIDF(h3k4me3)
h3k4me3 <- FindTopFeatures(h3k4me3, min.cutoff = '5')
h3k4me3 <- RunSVD(object = h3k4me3,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
pdf("h3k4me3_peaks.pdf")
DepthCor(h3k4me3)
h3k4me3 <- RunUMAP(object = h3k4me3,reduction = 'lsi',dims = c(2:15))
h3k4me3 <- FindNeighbors(object = h3k4me3,reduction = 'lsi',dims = c(2:15))
h3k4me3 <- FindClusters(object = h3k4me3,algorithm = 1,resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k4me3, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k4me3@active.ident))
label <- as.data.frame(label)
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("16cell",label[,1]),2] <- "morula"
label[grep("^T5",label[,1]),2] <- "blastocyst"
label[,3] <- h3k4me3@meta.data$seurat_clusters
h3k4me3@meta.data$group.ident <- as.factor(label[,2])
DimPlot(object = h3k4me3, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()

#################### h3k27ac ########################
pathToBams1 <- c('../2cell/h3k27ac/')
bamFiles1 <- paste(pathToBams1, inte, sep='')
pathToBams2 <- c('../4cell/h3k27ac/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
pathToBams3 <- c('../8cell/h3k27ac/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
pathToBams4 <- c('../../05_morula_inte/h3k27ac/')
bamFiles4 <- paste(pathToBams4, list.files(pathToBams4), sep='')
bamFiles4 <- bamFiles4[grep(".bam", bamFiles4)]
pathToBams5 <- c('../../06_blasto_inte/03_scChromHMM/02_chromHMM_sc/bam/h3k27ac/')
bamFiles5 <- paste(pathToBams5, list.files(pathToBams5), sep='')
pathToBams6 <- c('../zygote/h3k27ac/')
bamFiles6 <- paste(pathToBams6, list.files(pathToBams6), sep='')
bamFiles <- c(bamFiles1, bamFiles2, bamFiles3, bamFiles4, bamFiles5,bamFiles6)

regions <- "./mm10.10K.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, paired = T, project.name='h3k27ac_linked', min.regions = 0)
c <- as.matrix(cisTopicObject@count.matrix) 

h3k27ac <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k27ac',min.cells = 1)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = '5')
h3k27ac <- RunSVD(object = h3k27ac,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
pdf("h3k27ac_10kbbin.pdf")
DepthCor(h3k27ac)
h3k27ac <- RunUMAP(object = h3k27ac,reduction = 'lsi',dims = c(2:15))
h3k27ac <- FindNeighbors(object = h3k27ac,reduction = 'lsi',dims = c(2:15))
h3k27ac <- FindClusters(object = h3k27ac,algorithm = 1,resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k27ac, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k27ac@active.ident))
label <- as.data.frame(label)
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("16cell",label[,1]),2] <- "morula"
label[grep("^T5",label[,1]),2] <- "blastocyst"
label[,3] <- h3k27ac@meta.data$seurat_clusters
h3k27ac@meta.data$group.ident <- as.factor(label[,2])
DimPlot(object = h3k27ac, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()


#################### h3k36me3 ########################
pathToBams1 <- c('../2cell/h3k36me3/')
bamFiles1 <- paste(pathToBams1, inte, sep='')
pathToBams2 <- c('../4cell/h3k36me3/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
pathToBams3 <- c('../8cell/h3k36me3/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
pathToBams4 <- c('../../05_morula_inte/h3k36me3/')
bamFiles4 <- paste(pathToBams4, list.files(pathToBams4), sep='')
bamFiles4 <- bamFiles4[grep(".bam", bamFiles4)]
pathToBams5 <- c('../../06_blasto_inte/03_scChromHMM/02_chromHMM_sc/bam/h3k36me3/')
bamFiles5 <- paste(pathToBams5, list.files(pathToBams5), sep='')
pathToBams6 <- c('../zygote/h3k36me3/')
bamFiles6 <- paste(pathToBams6, list.files(pathToBams6), sep='')
bamFiles <- c(bamFiles1, bamFiles2, bamFiles3, bamFiles4, bamFiles5,bamFiles6)

regions <- "./mm10_Refseq.transcript.uniq2.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, paired = T, project.name='h3k36me3_linked', min.regions = 0)
c <- as.matrix(cisTopicObject@count.matrix) 

h3k36me3 <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k36me3',min.cells = 1)
h3k36me3 <- RunTFIDF(h3k36me3)
h3k36me3 <- FindTopFeatures(h3k36me3, min.cutoff = '5')
h3k36me3 <- RunSVD(object = h3k36me3,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
pdf("h3k36me3_genebody.pdf")
DepthCor(h3k36me3)
h3k36me3 <- RunUMAP(object = h3k36me3,reduction = 'lsi',dims = c(2:10))
h3k36me3 <- FindNeighbors(object = h3k36me3,reduction = 'lsi',dims = c(2:10))
h3k36me3 <- FindClusters(object = h3k36me3,algorithm = 1,resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k36me3, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k36me3@active.ident))
label <- as.data.frame(label)
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("16cell",label[,1]),2] <- "morula"
label[grep("^T5",label[,1]),2] <- "blastocyst"
label[,3] <- h3k36me3@meta.data$seurat_clusters
h3k36me3@meta.data$group.ident <- as.factor(label[,2])
DimPlot(object = h3k36me3, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()


#################### h3k27me3 ########################
pathToBams1 <- c('../2cell/h3k27me3/')
bamFiles1 <- paste(pathToBams1, inte, sep='')
pathToBams2 <- c('../4cell/h3k27me3/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
pathToBams3 <- c('../8cell/h3k27me3/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
pathToBams4 <- c('../../05_morula_inte/h3k27me3/')
bamFiles4 <- paste(pathToBams4, list.files(pathToBams4), sep='')
bamFiles4 <- bamFiles4[grep(".bam", bamFiles4)]
pathToBams5 <- c('../../06_blasto_inte/03_scChromHMM/02_chromHMM_sc/bam/h3k27me3/')
bamFiles5 <- paste(pathToBams5, list.files(pathToBams5), sep='')
pathToBams6 <- c('../zygote/h3k27me3/')
bamFiles6 <- paste(pathToBams6, list.files(pathToBams6), sep='')
bamFiles <- c(bamFiles1, bamFiles2, bamFiles3, bamFiles4, bamFiles5,bamFiles6)

regions <- "./mm10.10K.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, paired = T, project.name='h3k27me3_linked', min.regions = 0)
c <- as.matrix(cisTopicObject@count.matrix) 

h3k27me3 <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k27me3',min.cells = 1)
h3k27me3 <- RunTFIDF(h3k27me3)
h3k27me3 <- FindTopFeatures(h3k27me3, min.cutoff = '5')
h3k27me3 <- RunSVD(object = h3k27me3,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
pdf("h3k27me3_10kb.pdf")
DepthCor(h3k27me3)
h3k27me3 <- RunUMAP(object = h3k27me3,reduction = 'lsi',dims = c(2:20))
h3k27me3 <- FindNeighbors(object = h3k27me3,reduction = 'lsi',dims = c(2:20))
h3k27me3 <- FindClusters(object = h3k27me3,algorithm = 1,resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k27me3, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k27me3@active.ident))
label <- as.data.frame(label)
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("16cell",label[,1]),2] <- "morula"
label[grep("^T5",label[,1]),2] <- "blastocyst"
label[,3] <- h3k27me3@meta.data$seurat_clusters
h3k27me3@meta.data$group.ident <- as.factor(label[,2])
DimPlot(object = h3k27me3, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()


#################### h3k9me3 ########################
pathToBams1 <- c('../2cell/h3k9me3/')
bamFiles1 <- paste(pathToBams1, inte, sep='')
pathToBams2 <- c('../4cell/h3k9me3/')
bamFiles2 <- paste(pathToBams2, list.files(pathToBams2), sep='')
pathToBams3 <- c('../8cell/h3k9me3/')
bamFiles3 <- paste(pathToBams3, list.files(pathToBams3), sep='')
pathToBams4 <- c('../../05_morula_inte/h3k9me3/')
bamFiles4 <- paste(pathToBams4, list.files(pathToBams4), sep='')
bamFiles4 <- bamFiles4[grep(".bam", bamFiles4)]
pathToBams5 <- c('../../06_blasto_inte/03_scChromHMM/02_chromHMM_sc/bam/h3k9me3/')
bamFiles5 <- paste(pathToBams5, list.files(pathToBams5), sep='')
pathToBams6 <- c('../zygote/h3k9me3/')
bamFiles6 <- paste(pathToBams6, list.files(pathToBams6), sep='')
bamFiles <- c(bamFiles1, bamFiles2, bamFiles3, bamFiles4, bamFiles5,bamFiles6)

regions <- "./mm10.10K.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, paired = T, project.name='h3k9me3_linked', min.regions = 0)
c <- as.matrix(cisTopicObject@count.matrix) 

h3k9me3 <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k9me3',min.cells = 1)
h3k9me3 <- RunTFIDF(h3k9me3)
h3k9me3 <- FindTopFeatures(h3k9me3, min.cutoff = '5')
h3k9me3 <- RunSVD(object = h3k9me3,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
pdf("h3k9me3_10kb.pdf")
DepthCor(h3k9me3)
h3k9me3 <- RunUMAP(object = h3k9me3,reduction = 'lsi',dims = c(2:15))
h3k9me3 <- FindNeighbors(object = h3k9me3,reduction = 'lsi',dims = c(2:15))
h3k9me3 <- FindClusters(object = h3k9me3,algorithm = 1,resolution = 0.5,verbose = FALSE)
DimPlot(object = h3k9me3, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k9me3@active.ident))
label <- as.data.frame(label)
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("16cell",label[,1]),2] <- "morula"
label[grep("^T5",label[,1]),2] <- "blastocyst"
label[,3] <- h3k9me3@meta.data$seurat_clusters
h3k9me3@meta.data$group.ident <- as.factor(label[,2])
DimPlot(object = h3k9me3, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()

################################# WNN #############################
all <- CreateSeuratObject(
  counts = c,
  assay = 'h3k9me3',
  project = 'h3k9me3',
  min.cells = 1
)
DefaultAssay(all) <- "h3k9me3"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(object = all, assay = 'h3k9me3', reduction.key = 'h3k9me3lsi_', reduction.name = 'h3k9me3_lsi')
all <- RunUMAP(object = all, reduction = 'h3k9me3_lsi', dims = c(2:15),reduction.name = "umap.h3k9me3", reduction.key = "h3k9me3UMAP_")
all <- FindNeighbors(object = all, reduction = 'h3k9me3_lsi',dims = c(2:15))
all <- FindClusters(object = all,algorithm = 3, resolution = 0.5, verbose = FALSE)

all[["h3k4me1"]] <- h3k4me1@assays$peaks
all[["h3k4me3"]] <- h3k4me3@assays$peaks
all[["h3k27ac"]] <- h3k27ac@assays$peaks
all[["h3k27me3"]] <- h3k27me3@assays$peaks
all[["h3k36me3"]] <- h3k36me3@assays$peaks

DefaultAssay(all) <- "h3k4me1"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(object = all, assay = 'h3k4me1', reduction.key = 'h3k4me1lsi_', reduction.name = 'h3k4me1_lsi')
all <- RunUMAP(object = all, reduction = 'h3k4me1_lsi', dims = c(2:15), reduction.name = "umap.h3k4me1", reduction.key = "h3k4me1UMAP_")
all <- FindNeighbors(object = all, reduction = 'h3k4me1_lsi',dims = c(2:15))
all <- FindClusters(object = all,algorithm = 3, resolution = 0.5, verbose = FALSE)

DefaultAssay(all) <- "h3k4me3"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(object = all, assay = 'h3k4me3', reduction.key = 'h3k4me3lsi_', reduction.name = 'h3k4me3_lsi')
all <- RunUMAP(object = all, reduction = 'h3k4me3_lsi', dims = c(2:15), reduction.name = "umap.h3k4me3", reduction.key = "h3k4me3UMAP_")
all <- FindNeighbors(object = all, reduction = 'h3k4me3_lsi',dims = c(2:15))
all <- FindClusters(object = all,algorithm = 3, resolution = 0.5, verbose = FALSE)

DefaultAssay(all) <- "h3k27ac"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(object = all, assay = 'h3k27ac', reduction.key = 'h3k27aclsi_', reduction.name = 'h3k27ac_lsi')
all <- RunUMAP(object = all, reduction = 'h3k27ac_lsi', dims = c(2:15), reduction.name = "umap.h3k27ac", reduction.key = "h3k27acUMAP_")
all <- FindNeighbors(object = all, reduction = 'h3k27ac_lsi',dims = c(2:15))
all <- FindClusters(object = all,algorithm = 3, resolution = 0.5, verbose = FALSE)

DefaultAssay(all) <- "h3k27me3"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(object = all, assay = 'h3k27me3', reduction.key = 'h3k27me3lsi_', reduction.name = 'h3k27me3_lsi')
all <- RunUMAP(object = all, reduction = 'h3k27me3_lsi', dims = c(2:20), reduction.name = "umap.h3k27me3", reduction.key = "h3k27me3UMAP_")
all <- FindNeighbors(object = all, reduction = 'h3k27me3_lsi',dims = c(2:20))
all <- FindClusters(object = all,algorithm = 3, resolution = 0.5, verbose = FALSE)

DefaultAssay(all) <- "h3k36me3"
all <- RunTFIDF(all)
all <- FindTopFeatures(all, min.cutoff = 'q0')
all <- RunSVD(object = all, assay = 'h3k36me3', reduction.key = 'h3k36me3lsi_', reduction.name = 'h3k36me3_lsi')
all <- RunUMAP(object = all, reduction = 'h3k36me3_lsi', dims = c(2:15), reduction.name = "umap.h3k36me3", reduction.key = "h3k36me3UMAP_")
all <- FindNeighbors(object = all, reduction = 'h3k36me3_lsi',dims = c(2:15))
all <- FindClusters(object = all,algorithm = 3, resolution = 0.5, verbose = FALSE)


all <- FindMultiModalNeighbors(all, 
                                    reduction.list = list('h3k9me3_lsi', 'h3k4me1_lsi','h3k4me3_lsi', 'h3k27ac_lsi', 'h3k27me3_lsi','h3k36me3_lsi'), 
                                    dims.list = list(2:15, 2:15, 2:15, 2:15, 2:20, 2:15), 
                                    knn.range = 50)
all <- RunUMAP(all, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
all <- FindClusters(all, graph.name = "wsnn", algorithm = 3, resolution = 0.8, verbose = FALSE)
all$celltype <- all@meta.data$wsnn_res.0.8

pdf("WNN_cluster.pdf", width = 15, height = 10)
p1 <- DimPlot(all, reduction = "umap.h3k4me1", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k4me1")
p2 <- DimPlot(all, reduction = "umap.h3k4me3", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k4me3")
p3 <- DimPlot(all, reduction = "umap.h3k27ac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k27ac")
p4 <- DimPlot(all, reduction = "umap.h3k36me3", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k36me3")
p5 <- DimPlot(all, reduction = "umap.h3k27me3", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k27me3")
p6 <- DimPlot(all, reduction = "umap.h3k9me3", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k9me3")
p7 <- DimPlot(all, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2+p3+p4+p5+p6+p7+ plot_layout(ncol = 3) & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("WNN_stage.pdf", width = 15, height = 10)
p1 <- DimPlot(all, reduction = "umap.h3k4me1", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k4me1")
p2 <- DimPlot(all, reduction = "umap.h3k4me3", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k4me3")
p3 <- DimPlot(all, reduction = "umap.h3k27ac", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k27ac")
p4 <- DimPlot(all, reduction = "umap.h3k36me3", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k36me3")
p5 <- DimPlot(all, reduction = "umap.h3k27me3", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k27me3")
p6 <- DimPlot(all, reduction = "umap.h3k9me3", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("h3k9me3")
p7 <- DimPlot(all, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2+p3+p4+p5+p6+p7+ plot_layout(ncol = 3) & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("WNN_weight.pdf", width = 15, height = 10)
p1 <- FeaturePlot(object = all,reduction = "wnn.umap", features = "h3k4me1.weight",pt.size = 1.5,max.cutoff = 'q95')+ ggtitle("h3k4me1")
p2 <- FeaturePlot(object = all,reduction = "wnn.umap", features = "h3k4me3.weight",pt.size = 1.5,max.cutoff = 'q95')+ ggtitle("h3k4me3")
p3 <- FeaturePlot(object = all,reduction = "wnn.umap", features = "h3k27ac.weight",pt.size = 1.5,max.cutoff = 'q95')+ ggtitle("h3k27ac")
p4 <- FeaturePlot(object = all,reduction = "wnn.umap", features = "h3k36me3.weight",pt.size = 1.5,max.cutoff = 'q95')+ ggtitle("h3k36me3")
p5 <- FeaturePlot(object = all,reduction = "wnn.umap", features = "h3k27me3.weight",pt.size = 1.5,max.cutoff = 'q95')+ ggtitle("h3k27me3")
p6 <- FeaturePlot(object = all,reduction = "wnn.umap", features = "h3k9me3.weight",pt.size = 1.5,max.cutoff = 'q95')+ ggtitle("h3k9me3")
p7 <- DimPlot(all, reduction = "wnn.umap", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2+p3+p4+p5+p6+p7+ plot_layout(ncol = 3)  & theme(plot.title = element_text(hjust = 0.5))
dev.off()




