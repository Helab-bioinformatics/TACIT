
############################## Intergrate with RNA for 2-cell stage ###########################################
library(cisTopic)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(dplyr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(edgeR)
set.seed(1)

####################################### 2cell integration (show H3K4me3 as example) #############################
RNA_pseudo <- readRDS("./scRNA/RNA_pseudo.rds")
RNA_pseudo <- subset(RNA_pseudo,cells = colnames(RNA_pseudo[,grep("2cell", RNA_pseudo@meta.data$group.ident)]))
RNA_pseudo <- NormalizeData(RNA_pseudo)
RNA_pseudo <- FindVariableFeatures(RNA_pseudo, nfeatures = 10000)
RNA_pseudo[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo, pattern = "^MT-")
RNA_pseudo <- ScaleData(RNA_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_pseudo <- RunPCA(RNA_pseudo, npcs = 10)
RNA_pseudo <- RunUMAP(RNA_pseudo, dims = 1:10, reduction = 'pca',n.neighbors = 10L)
RNA_pseudo <- FindNeighbors(RNA_pseudo, dims = 1:10, reduction = "pca")
RNA_pseudo <- FindClusters(RNA_pseudo, resolution = 1)
pdf("./scRNA/RNA_pseudo_2cell.pdf")
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
dev.off()

label <- rownames(as.data.frame(RNA_pseudo@active.ident))
label <- as.data.frame(label)
label[,2] <- RNA_pseudo@meta.data$seurat_clusters

RNA_count <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_count <- RNA_count[which(rowSums(RNA_count) >0),]
toti <- read.csv("./totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("./pluripotency.markers.csv", sep = ",", header = F)
ZGA <- read.xlsx("/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/MajorZGA_genes.xlsx")
ZGA <- na.omit(ZGA$Gene)
RNA_ZGA <- RNA_count[which(rownames(RNA_count) %in% ZGA),]
RNA_toti <- RNA_count[which(rownames(RNA_count) %in% toti$V1),]
RNA_pluri <- RNA_count[which(rownames(RNA_count) %in% pluri$V1),]
RNA_ZGA <- colSums(RNA_ZGA)/620
RNA_toti <- colSums(RNA_toti)/215
RNA_pluri <- colSums(RNA_pluri)/38

label$ZGA <- RNA_ZGA
label$toti <- RNA_toti
label$pluri <- RNA_pluri
pdf("./scRNA/RNA_2cell_zga_average.pdf")
ggplot(label,mapping = aes(x = V2, y = pluri))+
  geom_violin(aes(fill =V2), trim = T) + 
  geom_boxplot(width = 0.05)+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0","#AC73EF"))+
  theme(legend.position = "none")+
  theme_bw()
dev.off()

RNA_pseudo@meta.data$zga.ident <- label$ZGA
RNA_pseudo@meta.data$zga.ident <- label$toti
RNA_pseudo@meta.data$zga.ident <- label$pluri

saveRDS(RNA_pseudo, "./scRNA/RNA_pseudo_2cell.rds")

############################################# H3K4me3 ####################################################
setwd("./h3k4me3")
h3k4me3 <- readRDS("./h3k4me3_allpeaks.rds")
h3k4me3 <- subset(h3k4me3,cells = colnames(h3k4me3[,grep("2cell", h3k4me3@meta.data$group.ident)]))
h3k4me3 <- RunTFIDF(h3k4me3)
h3k4me3 <- FindTopFeatures(h3k4me3, min.cutoff = 'q0')
h3k4me3 <- RunSVD(
  object = h3k4me3,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
h3k4me3 <- RunUMAP(object = h3k4me3, reduction = 'lsi', n.neighbors = 20L, dims = c(2:15))
h3k4me3 <- FindNeighbors(object = h3k4me3, reduction = 'lsi', dims = c(2:15))
h3k4me3 <- FindClusters(object = h3k4me3, algorithm = 3, resolution = 1, verbose = FALSE)
h3k4me3_count <- as.matrix(h3k4me3@assays[["peaks"]]@counts)
depth <- apply(h3k4me3_count, 2, sum)
h3k4me3$depth <- log10(depth)

pdf("./h3k4me3_2cell.pdf")
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

write.table(h3k4me3_pre1, "../h3k4me3/h3k4me3_2cell_pseudo_inte.txt", quote = F, sep = " ")
saveRDS(h3k4me3, "../h3k4me3/h3k4me3_2cell_inte.rds")
saveRDS(transfer.anchors, "../h3k4me3/h3k4me3_2cell_transfer_anchors.rds")

