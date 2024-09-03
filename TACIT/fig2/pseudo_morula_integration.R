

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

######################## 直接看single cell和pseudo-cells，能否分成ICM和TE potential，效果一般 #########################
RNA <- readRDS("./scRNA/all_RNA_no_outlier.rds")
RNA <- subset(RNA,cells = colnames(RNA[,grep("morula", colnames(RNA))]))
RNA <- NormalizeData(RNA)
RNA <- FindVariableFeatures(RNA, nfeatures = 10000)
RNA[["percent.mt"]] <- PercentageFeatureSet(RNA, pattern = "^MT-")
RNA <- ScaleData(RNA, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA <- RunPCA(RNA, npcs = 10)
RNA <- RunUMAP(RNA, dims = 1:10, reduction = 'pca',n.neighbors = 10L)
RNA <- FindNeighbors(RNA, dims = 1:10, reduction = "pca")
RNA <- FindClusters(RNA, resolution = 1)

pdf("./scRNA/RNA_morula_markers.pdf")
DimPlot(RNA, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA, group.by = "group.ident" ,label = TRUE, pt.size = 2)
FeaturePlot(object = RNA,features = c("Sox2", "Nanog", "Pou5f1", "Klf4", "Otx2", "Gata4"),pt.size = 1.5,max.cutoff = 'q95', ncol = 2)
FeaturePlot(object = RNA,features = c("Cdx2", "Eomes", "Gata2", "Id2", "Tbx2", "Elf5"),pt.size = 1.5,max.cutoff = 'q95', ncol = 2)
dev.off()

#######################################################################################################################
####################################### morula integration ######################################################
RNA_pseudo <- readRDS("./scRNA/RNA_pseudo.rds")
RNA_pseudo <- subset(RNA_pseudo,cells = colnames(RNA_pseudo[,grep("morula", colnames(RNA_pseudo))]))
RNA_pseudo <- NormalizeData(RNA_pseudo)
RNA_pseudo <- FindVariableFeatures(RNA_pseudo, nfeatures = 10000)
RNA_pseudo[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo, pattern = "^MT-")
RNA_pseudo <- ScaleData(RNA_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_pseudo <- RunPCA(RNA_pseudo, npcs = 10)
RNA_pseudo <- RunUMAP(RNA_pseudo, dims = 1:10, reduction = 'pca',n.neighbors = 10L)
RNA_pseudo <- FindNeighbors(RNA_pseudo, dims = 1:10, reduction = "pca")
RNA_pseudo <- FindClusters(RNA_pseudo, resolution = 1)
pdf("./scRNA/RNA_pseudo_morula.pdf")
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
FeaturePlot(object = RNA_pseudo,features = c("Sox2", "Nanog", "Pou5f1", "Klf4", "Otx2", "Gata4"),pt.size = 1.5,max.cutoff = 'q95', ncol = 2)
FeaturePlot(object = RNA_pseudo,features = c("Cdx2", "Eomes", "Gata2", "Id2", "Tbx2", "Elf5"),pt.size = 1.5,max.cutoff = 'q95', ncol = 2)
dev.off()

saveRDS(RNA_pseudo, "./scRNA/RNA_pseudo_morula.rds")
######single cells和pseudo-cells均没有很好分开ICM-和TE-potential,所以不另外合成pseudo-cells################


############ hcluster morula stage ################
RNA_counts <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_counts <- as.matrix(RNA_counts)
RNA_counts <- t(RNA_counts)
out.dist <- dist(RNA_counts, method = "euclidean")
out.hclust <- hclust(out.dist, method = "ward.D2")
out.id <- cutree(out.hclust, k=2)
#identify late 2-cell as a cluster
pdf("./scRNA/RNA_morula_pseudo.hcluster.pdf", width = 12)
plot(out.hclust)
dev.off()

############################################# H3K4me3 ####################################################
setwd("./h3k4me3")
h3k4me3 <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allpeaks.rds")
h3k4me3 <- subset(h3k4me3,cells = colnames(h3k4me3[,grep("morula", h3k4me3@meta.data$group.ident)]))
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

pdf("./h3k4me3_morula.pdf")
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

write.table(h3k4me3_pre1, "../h3k4me3/h3k4me3_morula_pseudo_inte.txt", quote = F, sep = " ")
saveRDS(h3k4me3, "../h3k4me3/h3k4me3_morula_inte.rds")
saveRDS(transfer.anchors, "../h3k4me3/h3k4me3_morula_transfer_anchors.rds")

######################################### H3K4me1 #############################################
############H3K4me1的4细胞数目少且reads数多，对应效果好，不需要构建pseudo-cells
setwd("../h3k4me1/")
set.seed(1)
pathToBams <- '/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/'
bamFiles <- paste(pathToBams, list.files(pathToBams), sep='')
bamFiles <- bamFiles[grep("morula", bamFiles)]
bamFiles <- bamFiles[grep("bai", invert = T, bamFiles)]

regions <- '/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/mm10.5K.windows.bed'
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles, regions, 
                                              min.cells = 0,#not filtering features
                                              allowMultiOverlap=T,#logical indicating if a read is allowed to be assigned to more than one feature (or meta-feature) if it is found to overlap with more than one feature (or meta-feature)
                                              project.name='H3K4me1',paired = T)
h3k4me1_counts <- as.matrix(cisTopicObject@count.matrix)

h3k4me1 <- CreateSeuratObject(counts = h3k4me1_counts,assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = 'q0')
h3k4me1 <- RunSVD(object = h3k4me1, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

DefaultAssay(h3k4me1) <- "peaks"
pdf("h3k4me1_morula_5kbbin.pdf")
DepthCor(h3k4me1)
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:10), n.neighbors = 10L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) 
label <- rownames(as.data.frame(h3k4me1@active.ident))
label <- as.data.frame(label)
label[grep("morula",label[,1]),2] <- "morula"
h3k4me1@meta.data$group.ident <- as.factor(label[,2])
DimPlot(h3k4me1, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off()

#geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me1/h3k4me1_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- read.csv("./h3k4me1_morula_allcells_tss15kb_ActivityMatrix.txt", sep = " ",check.names=F)
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
                                        k.anchor = 10, 
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
write.table(h3k4me1_pre1, "../h3k4me1/h3k4me1_morula_pseudo.txt", quote = F, sep = " ")

saveRDS(h3k4me1, "../h3k4me1/h3k4me1_morula_integration.rds")
saveRDS(transfer.anchors, "../h3k4me1/h3k4me1_morula_transfer_anchors.rds")

######################################### H3K27ac #######################################
setwd("../h3k27ac/")
set.seed(2)
h3k27ac_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/h3k27ac.5kb.best.countMatrix.txt", sep =" ") 
h3k27ac_counts <- h3k27ac_counts[,grep("morula",colnames(h3k27ac_counts))] 
h3k27ac <- CreateSeuratObject(counts = h3k27ac_counts,assay = 'peaks',project = 'h3k27ac',min.cells = 1)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = 'q0')
h3k27ac <- RunSVD(object = h3k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k27ac/h3k27ac_allStage_peaks_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,grep("morula", colnames(geneActivity))]

dat <- data.frame(col=colnames(geneActivity), min="min")
dat <- separate(dat, col, c("v1","v2", "v3"), sep = "-")
dat <- unite(dat,col, c("v1","v2", "v3"), sep = ".")
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
                                        k.anchor = 10, 
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
write.table(h3k27ac_pre1, "./h3k27ac_morula_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k27ac, "h3k27ac_morula_integration.rds")
saveRDS(transfer.anchors, "h3k27ac_morula_transfer_anchors.rds")


############################################### H3K36me3 #########################################
#rm(list = ls())
setwd("../h3k36me3/")
set.seed(3)
h3k36me3 <- readRDS("./h3k36me3_5kb.rds")
h3k36me3 <- subset(h3k36me3,cells = colnames(h3k36me3[,grep("morula", h3k36me3@meta.data$group.ident)]))
h3k36me3 <- RunTFIDF(h3k36me3)
h3k36me3 <- FindTopFeatures(h3k36me3, min.cutoff = 'q0')
h3k36me3 <- RunSVD(object = h3k36me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
h3k36me3 <- RunUMAP(object = h3k36me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
h3k36me3 <- FindNeighbors(object = h3k36me3,reduction = 'lsi',dims = c(2:15))
h3k36me3 <- FindClusters(object = h3k36me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
pdf("h3k36me3_morula_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k36me3/h3k36me3_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,grep("morula", colnames(geneActivity))]
geneActivity <- geneActivity[,colnames(h3k36me3)]
h3k36me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k36me3) <- "ACTIVITY"
h3k36me3 <- FindVariableFeatures(h3k36me3, nfeatures = 10000)
h3k36me3 <- NormalizeData(h3k36me3)
h3k36me3 <- ScaleData(h3k36me3)

#RNA <- readRDS("../RNA_morula.rds")
#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k36me3,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 10, 
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
write.table(h3k36me3_pre1, "./h3k36me3_morula_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k36me3, "h3k36me3_morula_integration.rds")
saveRDS(transfer.anchors, "h3k36me3_morula_transfer_anchors.rds")

########################################### H3K27me3 #############################################
##借助coTACIT数据
setwd("../cotacit_data/")
cotacit <- readRDS("./cotacit_h3k27ac.rds")
cotacit_morula <- subset(cotacit,cells = colnames(cotacit[,grep("morula", cotacit@meta.data$group.ident)]))
geneActivity <- read.csv("./cotacit_h3k27ac_allstage_tss5kb_ActivityMatrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(cotacit_morula)]
cotacit_morula[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(cotacit_morula) <- "ACTIVITY"
cotacit_morula <- FindVariableFeatures(cotacit_morula, nfeatures = 10000)
cotacit_morula <- NormalizeData(cotacit_morula)
cotacit_morula <- ScaleData(cotacit_morula)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = cotacit_morula,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 10, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = cotacit_morula[["lsi"]], dims = 2:ncol(cotacit_morula[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_pre1 <- ID.predict
length(unique(cotacit_pre1$predicted.id))
#cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]
write.table(cotacit_pre1, "./cotacit_morula_h3k27ac_pseudo.txt", quote = F, sep = " ")
saveRDS(cotacit_morula, "cotacit_morula_h3k27ac_integration.rds")
saveRDS(transfer.anchors, "cotacit_morula_h3k27ac_transfer_anchors.rds")


####################################### all #######################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_morula_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_morula_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_morula_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_morula_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_morula_h3k27ac_pseudo.txt", sep = " ")

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

write.table(h3k4me1_pre1, "./h3k4me1_morula_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k4me3_pre1, "./h3k4me3_morula_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k27ac_pre1, "./h3k27ac_morula_pseudo_select.txt", quote = F, sep = " ")
write.table(h3k36me3_pre1, "./h3k36me3_morula_pseudo_select.txt", quote = F, sep = " ")
write.table(cotacit_pre1, "./cotacit_morula_pseudo_select.txt", quote = F, sep = " ")
##由于想尽可能多地保留整合后的细胞，所以整合到同一个RNA细胞的多个细胞也分开并保留，因此在predicted.id后面还增加了一个编号，这一步在excel中操作(按照predicted.id和score值进行sort后编号)。
h3k4me1_pre1 <- read.csv("./h3k4me1_morula_pseudo_select.txt", sep = "\t")
h3k4me3_pre1 <- read.csv("./h3k4me3_morula_pseudo_select.txt", sep = "\t")
h3k27ac_pre1 <- read.csv("./h3k27ac_morula_pseudo_select.txt", sep = "\t")
h3k36me3_pre1 <- read.csv("./h3k36me3_morula_pseudo_select.txt", sep = "\t")
cotacit_pre1 <- read.csv("./cotacit_morula_pseudo_select.txt", sep = "\t")

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))
inte
h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1$predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1$predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1$predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1$predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$predicted.id %in% inte),]

h3k4me1.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/", h3k4me1_pre1$X, " ./morula/", h3k4me1_pre1$predicted.id, "_inte.bam")
write.table(h3k4me1.sh, "./h3k4me1/h3k4me1_morula_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", h3k4me3_pre1$X, " ./morula/", h3k4me3_pre1$predicted.id, "_inte.bam")
write.table(h3k4me3.sh, "./h3k4me3/h3k4me3_morula_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k27ac_pre1 <- separate(h3k27ac_pre1, X, c("x1", "x2", "x3", "x4", "x5"), sep = "\\.")
h3k27ac_pre1 <- unite(h3k27ac_pre1, id, c("x1", "x2", "x3"), sep = "-")
h3k27ac_pre1 <- unite(h3k27ac_pre1, X, c("id", "x4", "x5"), sep = ".")
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", h3k27ac_pre1$X, " ./morula/", h3k27ac_pre1$predicted.id, "_inte.bam")
write.table(h3k27ac.sh, "./h3k27ac/h3k27ac_morula_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", h3k36me3_pre1$X, " ./morula/", h3k36me3_pre1$predicted.id, "_inte.bam")
write.table(h3k36me3.sh, "./h3k36me3/h3k36me3_morula_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

cotacit_pre1$ID <- cotacit_pre1$X
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_all_h3k9me3.mm10_rmdup.bam")
cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2, " ./morula/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k27me3/cotacit_morula_linked_h3k27me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3, " ./morula/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k9me3/cotacit_morula_linked_h3k9me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(cotacit_pre1[,c("predicted.id","prediction.score.max", "ID2", "ID3")], "./cotacit_data/integration/cotacit_morula_integration.txt")


############################################# cotacit + tacit h3k37me3 #########################################

############ cotacit h3k27me3 ################
library(Signac)
h3k27me3 <- readRDS("./h3k27me3_all stage_peaks.rds")
cotacit_h3k27me3 <- readRDS("./cotacit_data/cotacit_h3k27me3_peaks.rds")

h3k27me3$dataset <- "tacit"
cotacit_h3k27me3$dataset <- "cotacit"

h3k27me3_morula <- subset(h3k27me3,cells = colnames(h3k27me3)[grep("E2.5", colnames(h3k27me3))])
cotacit_h3k27me3_morula <- subset(cotacit_h3k27me3,cells = colnames(cotacit_h3k27me3)[grep("morula", cotacit@meta.data$group.ident)])
merge <- merge(h3k27me3_morula, cotacit_h3k27me3_morula)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_morula.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k27me3_morula, cotacit_h3k27me3_morula),
  anchor.features = rownames(h3k27me3_morula),
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

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_integration_morula.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k27me3_morula,
  query = cotacit_h3k27me3_morula,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k27me3_morula@meta.data$ID <- colnames(h3k27me3_morula)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k27me3_morula$ID,
                                   weight.reduction = cotacit_h3k27me3_morula[["lsi"]], dims = 2:ncol(cotacit_h3k27me3_morula[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k27me3_pre1 <- ID.predict
length(unique(cotacit_h3k27me3_pre1$predicted.id))

write.table(cotacit_h3k27me3_pre1, "./cotacit_data/h3k27me3_morula_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k27me3_morula_cotacit_tacit_inte_transfer_anchors.rds")

############################################# cotacit + tacit h3k9me3 #########################################

############ cotacit h3k9me3 ################
cotacit_h3k9me3 <- readRDS("./cotacit_data/cotacit_h3k9me3_10kbbin.rds")
h3k9me3 <- readRDS("./cotacit_data/h3k9me3_10kbbin.rds")
h3k9me3$dataset <- "tacit"
cotacit_h3k9me3$dataset <- "cotacit"

h3k9me3_morula <- subset(h3k9me3,cells = colnames(h3k9me3)[grep("morula", h3k9me3@meta.data$group.ident)])
cotacit_h3k9me3_morula <- subset(cotacit_h3k9me3,cells = colnames(cotacit_h3k9me3)[grep("morula", cotacit_h3k9me3@meta.data$group.ident)])
merge <- merge(h3k9me3_morula, cotacit_h3k9me3_morula)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_morula.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k9me3_morula, cotacit_h3k9me3_morula),
  anchor.features = rownames(h3k9me3_morula),
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

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_integration_morula.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k9me3_morula,
  query = cotacit_h3k9me3_morula,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k9me3_morula@meta.data$ID <- colnames(h3k9me3_morula)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k9me3_morula$ID,
                                   weight.reduction = cotacit_h3k9me3_morula[["lsi"]], dims = 2:ncol(cotacit_h3k9me3_morula[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k9me3_pre1 <- ID.predict
length(unique(cotacit_h3k9me3_pre1$predicted.id))

write.table(cotacit_h3k9me3_pre1, "./cotacit_data/h3k9me3_morula_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k9me3_morula_cotacit_tacit_inte_transfer_anchors.rds")

cotacit_pre1$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1$ID3, rownames(cotacit_h3k9me3_pre1))]
k27me3.sh <- paste0("samtools merge ./h3k27me3_tacit+cotacit/morula/", cotacit_pre1$predicted.id, "_inte.bam ", 
                    "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2,
                    " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1$k27_tacit)
k9me3.sh <- paste0("samtools merge ./h3k9me3_tacit+cotacit/morula/", cotacit_pre1$predicted.id, "_inte.bam ", 
                   "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3,
                   " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1$k9_tacit)
write.table(k27me3.sh, "./cotacit_data/integration/h3k27me3_morula_merge.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./cotacit_data/integration/h3k9me3_morula_merge.sh", sep = "\t", quote = F, col.names = F, row.names = F)

write.table(cotacit_pre1, "./cotacit_data/cotacit_tacit_inte_morula.txt", quote = F, sep = " ")


################################## 为了用尽可能多的细胞构建chromHMM，基本保留所有整合的TACIT和cotacit细胞 ###########################################
###### 这里和pseudo_4cell_meta.R做的事情类似 #####
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_morula_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_morula_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_morula_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_morula_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_morula_h3k27ac_pseudo.txt", sep = " ")

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
write.table(h3k4me1.sh, "./08_chromHMM_morula/h3k4me1_morula_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k4me3 ####
h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", rownames(h3k4me3_pre1), 
                     " ./", h3k4me3_pre1$predicted.id, "/")
write.table(h3k4me3.sh, "./08_chromHMM_morula/h3k4me3_morula_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k27ac ####
h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k27ac_pre1$X <- rownames(h3k27ac_pre1)
h3k27ac_pre1 <- separate(h3k27ac_pre1, X, c("x1", "x2", "x3", "x4", "x5"), sep = "\\.")
h3k27ac_pre1 <- unite(h3k27ac_pre1, id, c("x1", "x2", "x3"), sep = "-")
h3k27ac_pre1 <- unite(h3k27ac_pre1, X, c("id", "x4", "x5"), sep = ".")
rownames(h3k27ac_pre1) <- h3k27ac_pre1$X
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", rownames(h3k27ac_pre1), 
                     " ./", h3k27ac_pre1$predicted.id, "/")
write.table(h3k27ac.sh, "./08_chromHMM_morula/h3k27ac_morula_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)

#### h3k36me3 ####
h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", rownames(h3k36me3_pre1), 
                      " ./", h3k36me3_pre1$predicted.id, "/")
write.table(h3k36me3.sh, "./08_chromHMM_morula/h3k36me3_morula_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)


#### cotacit ####
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_morula_h3k27ac_pseudo.txt", sep = " ")
cotacit_h3k27me3_pre1 <- read.csv("./cotacit_data/h3k27me3_morula_cotacit_tacit_inte.txt", sep = " ")
cotacit_h3k9me3_pre1 <- read.csv("./cotacit_data/h3k9me3_morula_cotacit_tacit_inte.txt", sep = " ")

cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1 $predicted.id %in% inte),]
cotacit_pre1$ID <- rownames(cotacit_pre1)
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_all_h3k9me3.mm10_rmdup.bam") 
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
write.table(k27me3.sh, "./08_chromHMM_morula/h3k27me3_morula_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./08_chromHMM_morula/h3k9me3_morula_meta_cp.sh", sep = "\t", quote = F, col.names = F, row.names = F)





