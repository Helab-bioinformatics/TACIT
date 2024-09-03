setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/01_allbin_umap/")

library(tidyr)
library(matrixStats)


#################################### enhancer ###################################
########### 对应上位置信息 ############
rm(list = ls())
library(tidyr)
filtered_mat_rou <- read.csv("./all_enhancer_strong_posterior_selected.txt", sep = "\t")
mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)

loci <- data.frame(loci=rownames(filtered_mat_rou), delete=filtered_mat_rou$deleted_rows)
loci$chr <- mm10$V1[match(loci$loci, mm10$V4)]
loci$start <- mm10$V2[match(loci$loci, mm10$V4)]
loci$end <- mm10$V3[match(loci$loci, mm10$V4)]
loci$end2 <- loci$start+200+(200*loci$delete)
loci <- unite(loci, loci2, c("chr", "start"), sep = "-", remove = F)
loci <- unite(loci, loci2, c("loci2", "end2"), sep = "-", remove = F)
rownames(filtered_mat_rou) <- loci$loci2
filtered_mat_rou$C2_0 <- filtered_mat_rou$C2_2
filtered_mat_rou <- filtered_mat_rou[,c("C2_0", colnames(filtered_mat_rou[,1:15]))]

################################### lsi umap ###########################
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

enhancer <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'enhancer',min.cells = 1)
enhancer <- RunTFIDF(enhancer)
enhancer <- FindTopFeatures(enhancer, min.cutoff = '5')
enhancer <- RunSVD(object = enhancer, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./enhancer_all_200bin_20240125.pdf")
DepthCor(enhancer)
enhancer <- RunUMAP(object = enhancer, reduction = 'lsi', n.neighbors = 5L, dims = c(1:7))#或可调试reduction
enhancer <- FindNeighbors(object = enhancer, reduction = 'lsi', dims = c(1:7))
enhancer <- FindClusters(object = enhancer, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
enhancer@meta.data[["ID"]] <- colnames(enhancer)
DimPlot(object = enhancer, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = enhancer, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = enhancer, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

da_peaks <- FindMarkers(
  object = enhancer,
  group.by = "orig.ident",
  ident.1 = c("C2"), 
  ident.2 = c("C8"),
  min.pct = 0.5,#设置越小，找的差异peaks越多，但也越不显著，均取前面部分即可。
  test.use = 'wilcox'
)

top <- rownames(da_peaks[which(da_peaks$p_val<0.1),])
top <- as.data.frame(top)
filtered_mat_rou <- filtered_mat_rou[top$top,]

################ top #######################
enhancer <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'enhancer',min.cells = 1)
enhancer <- RunTFIDF(enhancer)
enhancer <- FindTopFeatures(enhancer, min.cutoff = '5')
enhancer <- RunSVD(object = enhancer, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./enhancer_strong_200bin_top_20240125.pdf")
DepthCor(enhancer)
enhancer <- RunUMAP(object = enhancer, reduction = 'lsi', n.neighbors = 10L, dims = c(1:6))#或可调试reduction
enhancer <- FindNeighbors(object = enhancer, reduction = 'lsi', dims = c(1:6))
enhancer <- FindClusters(object = enhancer, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
enhancer@meta.data[["ID"]] <- colnames(enhancer)
DimPlot(object = enhancer, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = enhancer, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = enhancer, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()


top <- separate(top, top, c("v1", "v2", "v3"), sep = "-")
write.table(format(top,scientific=FALSE), "./enhancer_strong_C2vsC8_diff_20230622.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(enhancer, "./enhancer_200bp_top_20240125.rds")


#################################### promoter ###################################
########### 对应上位置信息 ############
rm(list = ls())
library(tidyr)
filtered_mat_rou <- read.csv("./all_promoter_posterior_selected.txt", sep = "\t")
mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)

loci <- data.frame(loci=rownames(filtered_mat_rou), delete=filtered_mat_rou$deleted_rows)
loci$chr <- mm10$V1[match(loci$loci, mm10$V4)]
loci$start <- mm10$V2[match(loci$loci, mm10$V4)]
loci$end <- mm10$V3[match(loci$loci, mm10$V4)]
loci$end2 <- loci$start+200+(200*loci$delete)
loci <- unite(loci, loci2, c("chr", "start"), sep = "-", remove = F)
loci <- unite(loci, loci2, c("loci2", "end2"), sep = "-", remove = F)
rownames(filtered_mat_rou) <- loci$loci2
filtered_mat_rou$C2_0 <- filtered_mat_rou$C2_2
filtered_mat_rou <- filtered_mat_rou[,c("C2_0", colnames(filtered_mat_rou[,1:15]))]

################################### lsi umap ###########################
promoter <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'promoter',min.cells = 1)
promoter <- RunTFIDF(promoter)
promoter <- FindTopFeatures(promoter, min.cutoff = '5')
promoter <- RunSVD(object = promoter, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./promoter_all_200bin_20240125.pdf")
DepthCor(promoter)
promoter <- RunUMAP(object = promoter, reduction = 'lsi', n.neighbors = 5L, dims = c(2:5))#或可调试reduction
promoter <- FindNeighbors(object = promoter, reduction = 'lsi', dims = c(2:5))
promoter <- FindClusters(object = promoter, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
promoter@meta.data[["ID"]] <- colnames(promoter)
DimPlot(object = promoter, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = promoter, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = promoter, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

da_peaks <- FindMarkers(
  object = promoter,
  group.by = "orig.ident",
  ident.1 = c("C2"), 
  ident.2 = c("C8"),
  min.pct = 0.5,#设置越小，找的差异peaks越多，但也越不显著，均取前面部分即可。
  test.use = 'wilcox'
)

top <- rownames(da_peaks[which(da_peaks$p_val<0.01),])
top <- as.data.frame(top)
filtered_mat_rou <- filtered_mat_rou[top$top,]

################ top #######################
promoter <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'promoter',min.cells = 1)
promoter <- RunTFIDF(promoter)
promoter <- FindTopFeatures(promoter, min.cutoff = '5')
promoter <- RunSVD(object = promoter, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./promoter_strong_200bin_top_20240125.pdf")
DepthCor(promoter)
promoter <- RunUMAP(object = promoter, reduction = 'lsi', n.neighbors = 5L, dims = c(1:5))#或可调试reduction
promoter <- FindNeighbors(object = promoter, reduction = 'lsi', dims = c(1:5))
promoter <- FindClusters(object = promoter, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
promoter@meta.data[["ID"]] <- colnames(promoter)
DimPlot(object = promoter, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = promoter, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = promoter, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

top <- separate(top, top, c("v1", "v2", "v3"), sep = "-")
write.table(format(top,scientific=FALSE), "./promoter_strong_C2vsC8_diff_20230622.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(promoter, "./promoter_200bp_top_20240125.rds")


############################################ transcription #########################################
rm(list = ls())
library(tidyr)
filtered_mat_rou <- read.csv("./all_transcription_posterior_selected.txt", sep = "\t")
mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)

loci <- data.frame(loci=rownames(filtered_mat_rou), delete=filtered_mat_rou$deleted_rows)
loci$chr <- mm10$V1[match(loci$loci, mm10$V4)]
loci$start <- mm10$V2[match(loci$loci, mm10$V4)]
loci$end <- mm10$V3[match(loci$loci, mm10$V4)]
loci$end2 <- loci$start+200+(200*loci$delete)
loci <- unite(loci, loci2, c("chr", "start"), sep = "-", remove = F)
loci <- unite(loci, loci2, c("loci2", "end2"), sep = "-", remove = F)
rownames(filtered_mat_rou) <- loci$loci2
filtered_mat_rou$C2_0 <- filtered_mat_rou$C2_2
filtered_mat_rou <- filtered_mat_rou[,c("C2_0", colnames(filtered_mat_rou[,1:15]))]

################################### lsi umap ###########################
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(tidyr)

transcription <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'transcription',min.cells = 1)
transcription <- RunTFIDF(transcription)
transcription <- FindTopFeatures(transcription, min.cutoff = '5')
transcription <- RunSVD(object = transcription, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./transcription_all_200bin_20240125.pdf")
DepthCor(transcription)
transcription <- RunUMAP(object = transcription, reduction = 'lsi', n.neighbors = 10L, dims = c(1:5))#或可调试reduction
transcription <- FindNeighbors(object = transcription, reduction = 'lsi', dims = c(1:5))
transcription <- FindClusters(object = transcription, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
transcription@meta.data[["ID"]] <- colnames(transcription)
DimPlot(object = transcription, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = transcription, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = transcription, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

da_peaks <- FindMarkers(
  object = transcription,
  group.by = "orig.ident",
  ident.1 = c("C2"), 
  ident.2 = c("C8"),
  min.pct = 0.5,#设置越小，找的差异peaks越多，但也越不显著，均取前面部分即可。
  test.use = 'wilcox'
)

top <- rownames(da_peaks[which(da_peaks$p_val<0.01),])
top <- as.data.frame(top)
filtered_mat_rou <- filtered_mat_rou[top$top,]
filtered_mat_rou <- filtered_mat_rou[1:40000,]
################ top #######################
counts <- filtered_mat_rou[,2:15]
transcription <- CreateSeuratObject(counts = counts, assay = 'peaks', project = 'transcription',min.cells = 1)
transcription <- RunTFIDF(transcription)
transcription <- FindTopFeatures(transcription, min.cutoff = '5')
transcription <- RunSVD(object = transcription, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./transcription_strong_200bin_top_20240125.pdf")
DepthCor(transcription)
transcription <- RunUMAP(object = transcription, reduction = 'lsi', n.neighbors = 6L, dims = c(1:5))#或可调试reduction
transcription <- FindNeighbors(object = transcription, reduction = 'lsi', dims = c(1:5))
transcription <- FindClusters(object = transcription, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
transcription@meta.data[["ID"]] <- colnames(transcription)
DimPlot(object = transcription, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = transcription, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = transcription, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()


top <- separate(top, top, c("v1", "v2", "v3"), sep = "-")
write.table(format(top,scientific=FALSE), "./transcription_strong_C2vsC8_diff_20230622.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(transcription, "./transcription_200bp_top_20240125.rds")


################################### h3k27me3 ###################################
########### 对应上位置信息 ############
rm(list = ls())
filtered_mat_rou <- read.csv("./all_heterochromatin_h3k27me3_posterior_selected.txt", sep = "\t")
mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)

loci <- data.frame(loci=rownames(filtered_mat_rou), delete=filtered_mat_rou$deleted_rows)
loci$chr <- mm10$V1[match(loci$loci, mm10$V4)]
loci$start <- mm10$V2[match(loci$loci, mm10$V4)]
loci$end <- mm10$V3[match(loci$loci, mm10$V4)]
loci$end2 <- loci$start+200+(200*loci$delete)
loci <- unite(loci, loci2, c("chr", "start"), sep = "-", remove = F)
loci <- unite(loci, loci2, c("loci2", "end2"), sep = "-", remove = F)
rownames(filtered_mat_rou) <- loci$loci2
filtered_mat_rou$C2_0 <- filtered_mat_rou$C2_2
filtered_mat_rou <- filtered_mat_rou[,c("C2_0", colnames(filtered_mat_rou[,1:15]))]

################################### lsi umap ###########################
h3k27me3 <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'h3k27me3',min.cells = 1)
h3k27me3 <- RunTFIDF(h3k27me3)
h3k27me3 <- FindTopFeatures(h3k27me3, min.cutoff = '5')
h3k27me3 <- RunSVD(object = h3k27me3, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./h3k27me3_all_200bin_20240125.pdf")
DepthCor(h3k27me3)
h3k27me3 <- RunUMAP(object = h3k27me3, reduction = 'lsi', n.neighbors = 10L, dims = c(1:5))#或可调试reduction
h3k27me3 <- FindNeighbors(object = h3k27me3, reduction = 'lsi', dims = c(1:5))
h3k27me3 <- FindClusters(object = h3k27me3, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
h3k27me3@meta.data[["ID"]] <- colnames(h3k27me3)
DimPlot(object = h3k27me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = h3k27me3, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = h3k27me3, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

da_peaks <- FindMarkers(
  object = h3k27me3,
  group.by = "orig.ident",
  ident.1 = c("C2"), 
  ident.2 = c("C8"),
  min.pct = 0.5,#设置越小，找的差异peaks越多，但也越不显著，均取前面部分即可。
  test.use = 'wilcox'
)

top <- rownames(da_peaks[which(da_peaks$p_val<0.01),])
top <- as.data.frame(top)
filtered_mat_rou <- filtered_mat_rou[top$top,]

################ top #######################
h3k27me3 <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'h3k27me3',min.cells = 1)
h3k27me3 <- RunTFIDF(h3k27me3)
h3k27me3 <- FindTopFeatures(h3k27me3, min.cutoff = '5')
h3k27me3 <- RunSVD(object = h3k27me3, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./h3k27me3_200bin_top_20240125.pdf")
DepthCor(h3k27me3)
h3k27me3 <- RunUMAP(object = h3k27me3, reduction = 'lsi', n.neighbors = 5L, dims = c(1:5))#或可调试reduction
h3k27me3 <- FindNeighbors(object = h3k27me3, reduction = 'lsi', dims = c(1:5))
h3k27me3 <- FindClusters(object = h3k27me3, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
h3k27me3@meta.data[["ID"]] <- colnames(h3k27me3)
DimPlot(object = h3k27me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = h3k27me3, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = h3k27me3, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()


top <- separate(top, top, c("v1", "v2", "v3"), sep = "-")
write.table(format(top,scientific=FALSE), "./h3k27me3_C2vsC8_diff_20230622.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(h3k27me3, "./h3k27me3_200bp_top_20240125.rds")



################################### h3k9me3 ###################################
########### 对应上位置信息 ############
rm(list = ls())
library(tidyr)
filtered_mat_rou <- read.csv("./all_heterochromatin_h3k9me3_posterior_selected.txt", sep = "\t")
mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)

loci <- data.frame(loci=rownames(filtered_mat_rou), delete=filtered_mat_rou$deleted_rows)
loci$chr <- mm10$V1[match(loci$loci, mm10$V4)]
loci$start <- mm10$V2[match(loci$loci, mm10$V4)]
loci$end <- mm10$V3[match(loci$loci, mm10$V4)]
loci$end2 <- loci$start+200+(200*loci$delete)
loci <- unite(loci, loci2, c("chr", "start"), sep = "-", remove = F)
loci <- unite(loci, loci2, c("loci2", "end2"), sep = "-", remove = F)
rownames(filtered_mat_rou) <- loci$loci2
filtered_mat_rou$C2_0 <- filtered_mat_rou$C2_2
filtered_mat_rou <- filtered_mat_rou[,c("C2_0", colnames(filtered_mat_rou[,1:15]))]

################################### lsi umap ###########################
h3k9me3 <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'h3k9me3',min.cells = 1)
h3k9me3 <- RunTFIDF(h3k9me3)
h3k9me3 <- FindTopFeatures(h3k9me3, min.cutoff = '5')
h3k9me3 <- RunSVD(object = h3k9me3, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./h3k9me3_all_200bin_20240125.pdf")
DepthCor(h3k9me3)
h3k9me3 <- RunUMAP(object = h3k9me3, reduction = 'lsi', n.neighbors = 10L, dims = c(1:5))#或可调试reduction
h3k9me3 <- FindNeighbors(object = h3k9me3, reduction = 'lsi', dims = c(1:5))
h3k9me3 <- FindClusters(object = h3k9me3, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
h3k9me3@meta.data[["ID"]] <- colnames(h3k9me3)
DimPlot(object = h3k9me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = h3k9me3, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = h3k9me3, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

da_peaks <- FindMarkers(
  object = h3k9me3,
  group.by = "orig.ident",
  ident.1 = c("C2"), 
  ident.2 = c("C8"),
  min.pct = 0.5,#设置越小，找的差异peaks越多，但也越不显著，均取前面部分即可。
  test.use = 'wilcox'
)

write.table(da_peaks, "./h3k9me3_C2vsC8_diff_all_20240126.txt", sep = "\t", quote = F, col.names = F, row.names = F)


top <- rownames(da_peaks[which(da_peaks$p_val<0.1),])
top <- as.data.frame(top)
filtered_mat_rou <- filtered_mat_rou[top$top,]
filtered_mat_rou <- filtered_mat_rou[1:3000,]
################ top #######################
h3k9me3 <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'h3k9me3',min.cells = 1)
h3k9me3 <- RunTFIDF(h3k9me3)
h3k9me3 <- FindTopFeatures(h3k9me3, min.cutoff = '5')
h3k9me3 <- RunSVD(object = h3k9me3, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./h3k9me3_200bin_top_20240125.pdf")
DepthCor(h3k9me3)
h3k9me3 <- RunUMAP(object = h3k9me3, reduction = 'lsi', n.neighbors = 9L, dims = c(1:10))#或可调试reduction
h3k9me3 <- FindNeighbors(object = h3k9me3, reduction = 'lsi', dims = c(1:10))
h3k9me3 <- FindClusters(object = h3k9me3, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
h3k9me3@meta.data[["ID"]] <- colnames(h3k9me3)
DimPlot(object = h3k9me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = h3k9me3, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = h3k9me3, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

top <- separate(top, top, c("v1", "v2", "v3"), sep = "-")
write.table(format(top,scientific=FALSE), "./h3k9me3_C2vsC8_diff_20230622.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(h3k9me3, "./h3k9me3_200bp_top_20240125.rds")



################################### heterochromatin ###################################
########### 对应上位置信息 ############
rm(list = ls())
library(tidyr)
filtered_mat_rou <- read.csv("./all_heterochromatin_both_posterior_selected.txt", sep = "\t")
mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)

loci <- data.frame(loci=rownames(filtered_mat_rou), delete=filtered_mat_rou$deleted_rows)
loci$chr <- mm10$V1[match(loci$loci, mm10$V4)]
loci$start <- mm10$V2[match(loci$loci, mm10$V4)]
loci$end <- mm10$V3[match(loci$loci, mm10$V4)]
loci$end2 <- loci$start+200+(200*loci$delete)
loci <- unite(loci, loci2, c("chr", "start"), sep = "-", remove = F)
loci <- unite(loci, loci2, c("loci2", "end2"), sep = "-", remove = F)
rownames(filtered_mat_rou) <- loci$loci2
filtered_mat_rou$C2_0 <- filtered_mat_rou$C2_2
filtered_mat_rou$C2_1 <- filtered_mat_rou$C2_2
filtered_mat_rou <- filtered_mat_rou[,c("C2_0", "C2_1", colnames(filtered_mat_rou[,1:14]))]

################################### lsi umap ###########################
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)
library(tidyr)

heterochromatin <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'heterochromatin',min.cells = 1)
heterochromatin <- RunTFIDF(heterochromatin)
heterochromatin <- FindTopFeatures(heterochromatin, min.cutoff = '5')
heterochromatin <- RunSVD(object = heterochromatin, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./heterochromatin_all_200bin_20240125.pdf")
DepthCor(heterochromatin)
heterochromatin <- RunUMAP(object = heterochromatin, reduction = 'lsi', n.neighbors = 5L, dims = c(1:5))#或可调试reduction
heterochromatin <- FindNeighbors(object = heterochromatin, reduction = 'lsi', dims = c(1:5))
heterochromatin <- FindClusters(object = heterochromatin, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
heterochromatin@meta.data[["ID"]] <- colnames(heterochromatin)
DimPlot(object = heterochromatin, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = heterochromatin, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = heterochromatin, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()

da_peaks <- FindMarkers(
  object = heterochromatin,
  group.by = "orig.ident",
  ident.1 = c("C4", "C2"), 
  ident.2 = c("C8"),
  min.pct = 0.05,#设置越小，找的差异peaks越多，但也越不显著，均取前面部分即可。
  test.use = 'wilcox'
)

top <- rownames(da_peaks[1:50000,])
top <- as.data.frame(top)
filtered_mat_rou <- filtered_mat_rou[top$top,]

################ top #######################
heterochromatin <- CreateSeuratObject(counts = filtered_mat_rou[,1:ncol(filtered_mat_rou)-1], assay = 'peaks', project = 'heterochromatin',min.cells = 1)
heterochromatin <- RunTFIDF(heterochromatin)
heterochromatin <- FindTopFeatures(heterochromatin, min.cutoff = '5')
heterochromatin <- RunSVD(object = heterochromatin, assay = 'peaks', reduction.key = 'LSI_',reduction.name = 'lsi')

pdf("./heterochromatin_strong_200bin_top_20240125.pdf")
DepthCor(heterochromatin)
heterochromatin <- RunUMAP(object = heterochromatin, reduction = 'lsi', n.neighbors = 10L, dims = c(1:5))#或可调试reduction
heterochromatin <- FindNeighbors(object = heterochromatin, reduction = 'lsi', dims = c(1:5))
heterochromatin <- FindClusters(object = heterochromatin, algorithm = 1, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
heterochromatin@meta.data[["ID"]] <- colnames(heterochromatin)
DimPlot(object = heterochromatin, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = heterochromatin, group.by = "orig.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = heterochromatin, group.by = "ID" ,label = TRUE, pt.size = 2) 
dev.off()


top <- separate(top, top, c("v1", "v2", "v3"), sep = "-")
write.table(format(top,scientific=FALSE), "./heterochromatin_C2vsC8_diff_20230622.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(heterochromatin, "./heterochromatin_200bp_top_20240125.rds")



