

#######################################################################################################################
####################################### pseudo-cell integration #######################################################
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
set.seed(1)

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

############################################# H3K4me3 (as an example) ####################################################
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

geneActivity <- read.csv("./h3k4me3_allstage_ActivitygeneActivityrix.txt", sep = " ",check.names=F)
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

####################################### all  (only retain cells with prediction.score.max >0.2 across all modalities) #######################################################
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



