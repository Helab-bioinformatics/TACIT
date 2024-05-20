############################ stage umap for each histone markers based on peaks #####################################
rm(list = ls())  
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F) 
library(cisTopic)
library(ggpubr)
library(ggsignif)
library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(tidyr)
library(stringr)
library(Matrix)
library(igraph)


########################################### h3k4me1 (as an example) #######################################
setwd("./all_h3k4me1/")
### zygote ###
pathToBams1 <- c('./allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|bai",invert = T, bamFiles1)]
bamFiles1 <- bamFiles1[grep("zygote", bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_zygote-p1e-5_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
h3k4me1 <- CreateSeuratObject(counts = c, assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:20),n.neighbors = 10L)
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:20))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.6,verbose = FALSE)
pdf("./allBam/03_stage_umap/peaks/h3k4me1_zygote_peaks.pdf")
DepthCor(h3k4me1)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

### 2cell ###
pathToBams1 <- c('./allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|bai",invert = T, bamFiles1)]
bamFiles1 <- bamFiles1[grep("2cell", bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_2cell-p1e-2_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
h3k4me1 <- CreateSeuratObject(counts = c, assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:10),n.neighbors = 12L)
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.6,verbose = FALSE)
pdf("./allBam/03_stage_umap/peaks/h3k4me1_2cell_peaks.pdf")
DepthCor(h3k4me1)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

### 4cell ###
pathToBams1 <- c('./allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|bai",invert = T, bamFiles1)]
bamFiles1 <- bamFiles1[grep("4cell", bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_4cell-0.05_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
h3k4me1 <- CreateSeuratObject(counts = c, assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:10),n.neighbors = 12L)
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.6,verbose = FALSE)
pdf("./allBam/03_stage_umap/peaks/h3k4me1_4cell_peaks.pdf")
DepthCor(h3k4me1)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

### 8cell ###
pathToBams1 <- c('./allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|bai",invert = T, bamFiles1)]
bamFiles1 <- bamFiles1[grep("8cell", bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_8cell-0.05_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
h3k4me1 <- CreateSeuratObject(counts = c, assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(1:10),n.neighbors = 12L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(1:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.6,verbose = FALSE)#或可调试reduction和algorithm
pdf("./allBam/03_stage_umap/peaks/h3k4me1_8cell_peaks.pdf")
DepthCor(h3k4me1)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

### morula ###
pathToBams1 <- c('./allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|bai",invert = T, bamFiles1)]
bamFiles1 <- bamFiles1[grep("morula", bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_morula-0.05_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
h3k4me1 <- CreateSeuratObject(counts = c, assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(3:30),n.neighbors = 12L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(3:30))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.6,verbose = FALSE)#或可调试reduction和algorithm
pdf("./allBam/03_stage_umap/peaks/h3k4me1_morula_peaks.pdf")
DepthCor(h3k4me1)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

### blastocyst ###
pathToBams1 <- c('./allBam/Min/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|bai",invert = T, bamFiles1)]
bamFiles1 <- bamFiles1[grep("blastocyst", bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_blastocyst-0.05_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 
h3k4me1 <- CreateSeuratObject(counts = c, assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(object = h3k4me1,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:30),n.neighbors = 20L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:30))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 1,resolution = 0.6,verbose = FALSE)#或可调试reduction和algorithm
pdf("./allBam/03_stage_umap/peaks/h3k4me1_blastocyst_peaks.pdf")
DepthCor(h3k4me1)
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

