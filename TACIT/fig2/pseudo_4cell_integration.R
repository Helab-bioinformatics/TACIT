
############################## Intergrate with RNA for 4-cell stage ###########################################
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


#####用所有peaks构建的geneActivityMatrix效果比5kb bin的好
#####把所有stage一起与RNA整合后，发现有很多细胞都match错了stage，所以尝试分stage来match，分zygote_2cell, 4cell_8cell, morula_blastocyst分别与RNA整合。
#####换用新的RNA-seq数据后，RNA细胞太多，导致很难有多个marker同时map到同一个细胞的情况，所以需要先对RNA-seq数据构建pseudo-cells

############################################# RNA ########################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/")

RNA <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/all_RNA_no_outlier.rds")
RNA <- subset(RNA,cells = colnames(RNA[,grep("4cell", RNA@meta.data$group.ident)]))
RNA <- NormalizeData(RNA)
RNA <- FindVariableFeatures(RNA, nfeatures = 5000)
RNA <- ScaleData(RNA, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA <- RunPCA(RNA, npcs = 10)
RNA <- RunUMAP(RNA, dims = 1:10, reduction = 'pca',n.neighbors = 10L)
RNA <- FindNeighbors(RNA, dims = 1:10, reduction = "pca")
RNA <- FindClusters(RNA, resolution = 1)
pdf("./4cell_RNA.pdf")
DimPlot(RNA, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(RNA, reduction = "umap", group.by = "group.ident",label = TRUE, pt.size = 2)
FeaturePlot(RNA, reduction = "umap", features = 'nFeature_RNA')
dev.off()

########### monocle for defining pseudo cells
RNA <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/all_RNA_no_outlier.rds")

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
#umap降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('cds.umap')
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(RNA, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="stage", cell_size = 0.8) + ggtitle('int.umap')
pdf("RNA_trajectory.pdf")
p1
p2
## 识别轨迹
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

#pesudo-cells, merge every 5 cells, for each stage one by one
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

saveRDS(RNA_pseudo, "../scRNA/RNA_pseudo.rds")
write.table(RNA_pseudo, "../scRNA/RNA_speudo_count.txt", quote = F, sep = "\t")

#######################################################################################################################
####################################### 4cell integration ######################################################
RNA_pseudo <- subset(RNA_pseudo,cells = colnames(RNA_pseudo[,grep("4cell", RNA_pseudo@meta.data$group.ident)]))
RNA_pseudo <- NormalizeData(RNA_pseudo)
RNA_pseudo <- FindVariableFeatures(RNA_pseudo, nfeatures = 10000)
RNA_pseudo[["percent.mt"]] <- PercentageFeatureSet(RNA_pseudo, pattern = "^MT-")
RNA_pseudo <- ScaleData(RNA_pseudo, vars.to.regress = c("nFeature_RNA", "percent.mito"))
RNA_pseudo <- RunPCA(RNA_pseudo, npcs = 10)
RNA_pseudo <- RunUMAP(RNA_pseudo, dims = 1:10, reduction = 'pca',n.neighbors = 5L)
RNA_pseudo <- FindNeighbors(RNA_pseudo, dims = 1:10, reduction = "pca")
RNA_pseudo <- FindClusters(RNA_pseudo, resolution = 1)

pdf("./RNA_pseudo_4cell.pdf")
DimPlot(RNA_pseudo, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = RNA_pseudo, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
dev.off()

label <- rownames(as.data.frame(RNA_pseudo@active.ident))
label <- as.data.frame(label)
label[,2] <- RNA_pseudo@meta.data$seurat_clusters

RNA_count <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_count <- RNA_count[which(rowSums(RNA_count) >0),]
toti <- read.csv("../totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("../pluripotency.markers.csv", sep = ",", header = F)
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
pdf("./RNA_4cell_zga_average.pdf")
ggplot(label,mapping = aes(x = V2, y = pluri))+
  geom_violin(aes(fill =V2), trim = T) + 
  geom_boxplot(width = 0.05)+
  scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0","#AC73EF"))+
  theme(legend.position = "none")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())
dev.off()

RNA_pseudo@meta.data$zga.ident <- label$ZGA
RNA_pseudo@meta.data$zga.ident <- label$toti
RNA_pseudo@meta.data$zga.ident <- label$pluri

saveRDS(RNA_pseudo, "./RNA_pseudo_4cell.rds")

############ hcluster 4cell stage ################
RNA_counts <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_counts <- as.matrix(RNA_counts)
RNA_counts <- t(RNA_counts)
toti <- read.csv("../totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("../pluripotency.markers.csv", sep = ",", header = F)

RNA_counts <- RNA_counts[,which(colnames(RNA_counts) %in% select)]
out.dist <- dist(RNA_counts, method = "manhattan")
out.hclust <- hclust(out.dist, method = "ward.D2")
out.id <- cutree(out.hclust, k=2)

pdf("./RNA_4cell_pseudo.hcluster.pdf")
plot(out.hclust)
dev.off()


############################################# H3K4me3 ####################################################

setwd("../h3k4me3")
h3k4me3 <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allpeaks.rds")
h3k4me3 <- subset(h3k4me3,cells = colnames(h3k4me3[,grep("4cell", h3k4me3@meta.data$group.ident)]))
h3k4me3 <- RunTFIDF(h3k4me3)
h3k4me3 <- FindTopFeatures(h3k4me3, min.cutoff = 'q0')
h3k4me3 <- RunSVD(
  object = h3k4me3,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
h3k4me3 <- RunUMAP(object = h3k4me3, reduction = 'lsi', dims = c(2:15))#或可调试reduction
h3k4me3 <- FindNeighbors(object = h3k4me3, reduction = 'lsi', dims = c(2:15))
h3k4me3 <- FindClusters(object = h3k4me3, algorithm = 3, resolution = 1, verbose = FALSE)#或可调试reduction和algorithm
h3k4me3_count <- as.matrix(h3k4me3@assays[["peaks"]]@counts)
depth <- apply(h3k4me3_count, 2, sum)
h3k4me3$depth <- log10(depth)

pdf("./h3k4me3_4cell.pdf")
DepthCor(h3k4me3)
DimPlot(object = h3k4me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(h3k4me3, reduction = "umap", group.by = "group.ident",label = TRUE, pt.size = 2)
FeaturePlot(h3k4me3, reduction = "umap", features = 'depth')
dev.off()

# geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allstage_ActivitygeneActivityrix.txt", sep = " ",check.names=F)
# geneActivity <- geneActivity[,grep("2cell|8cell|blasto", colnames(geneActivity))]
# h3k4me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
# DefaultAssay(h3k4me3) <- "ACTIVITY"
# h3k4me3 <- FindVariableFeatures(h3k4me3, nfeatures = 10000)
# h3k4me3 <- NormalizeData(h3k4me3)
# h3k4me3 <- ScaleData(h3k4me3)

###尝试promoter区域建立gene matrix
pathToBams1 <- c('/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bai|ArchR|cluster|bed",invert = T, bamFiles1)]
#regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl2kb_all+.bed"
regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl5kb_all+.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='h3k4me3', min.regions = 0)
h3k4me3_pr <- as.matrix(cisTopicObject@count.matrix)
##use CreatGeneActivityMatrix.R to generate gene activity matrix


geneActivity <- read.csv("./h3k4me3_allstage_tss2kb_ActivityMatrix.txt", sep = " ",check.names=F)
#geneActivity <- read.csv("./h3k4me3_allstage_tss5kb_ActivityMatrix.txt", sep = " ",check.names=F)
#geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allstage_ActivitygeneActivityrix.txt", sep = " ",check.names=F)
#geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me3/h3k4me3_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
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
                                        #                                      features = top500$gene,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 5)
#Annotate scATAC-seq cells via label transfer

#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = h3k4me3[["lsi"]], dims = 2:ncol(h3k4me3[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k4me3_pre1 <- ID.predict
length(unique(h3k4me3_pre1$predicted.id))
#h3k4me3_pre1 <- h3k4me3_pre1[which(h3k4me3_pre1$prediction.score.max >0.3),]

write.table(h3k4me3_pre1, "../h3k4me3/h3k4me3_4cell_pseudo_inte.txt", quote = F, sep = " ")
saveRDS(h3k4me3, "../h3k4me3/h3k4me3_4cell_inte.rds")
saveRDS(transfer.anchors, "../h3k4me3/h3k4me3_4cell_transfer_anchors.rds")

######################################### H3K4me1 #############################################
############H3K4me1的4细胞数目少且reads数多，对应效果好，不需要构建pseudo-cells
setwd("../h3k4me1/")
set.seed(1)
h3k4me1_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/h3k4me1.5kb.countMatrix.HX.txt", sep ="\t") 
rownames(h3k4me1_counts) <- h3k4me1_counts[,1]
h3k4me1_counts <- h3k4me1_counts[,-1]
h3k4me1_counts <- h3k4me1_counts[,grep("4cell",colnames(h3k4me1_counts))] 
h3k4me1 <- CreateSeuratObject(counts = h3k4me1_counts,assay = 'peaks',project = 'h3k4me1',min.cells = 1)
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = 'q0')
h3k4me1 <- RunSVD(object = h3k4me1, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

DefaultAssay(h3k4me1) <- "peaks"
pdf("h3k4me1_4cell_5kbbin.pdf")
DepthCor(h3k4me1)
h3k4me1 <- RunUMAP(object = h3k4me1,reduction = 'lsi',dims = c(2:10), n.neighbors = 10L)#或可调试reduction
h3k4me1 <- FindNeighbors(object = h3k4me1,reduction = 'lsi',dims = c(2:10))
h3k4me1 <- FindClusters(object = h3k4me1,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) 
label <- rownames(as.data.frame(h3k4me1@active.ident))
label <- as.data.frame(label)
label[grep("4cell",label[,1]),2] <- "4cell"
h3k4me1@meta.data$group.ident <- as.factor(label[,2])
DimPlot(h3k4me1, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off()

#geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me1/h3k4me1_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k4me1/h3k4me1_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(h3k4me1)]
h3k4me1[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k4me1) <- "ACTIVITY"
h3k4me1 <- FindVariableFeatures(h3k4me1, nfeatures = 10000)
h3k4me1 <- NormalizeData(h3k4me1)
h3k4me1 <- ScaleData(h3k4me1)

###尝试promoter区域建立gene matrix
# pathToBams1 <- c('/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/')
# bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
# bamFiles1 <- bamFiles1[grep("bai|ArchR|cluster|bed",invert = T, bamFiles1)]
# bamFiles1 <- bamFiles1[grep("4cell",bamFiles1)]
# #regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl2kb_all+.bed"
# regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl5kb_all+.bed"
# regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/05_coverage/agg/h3k4me1_4cell-p1e-2.80000.peaks.bed"
# cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='h3k4me3', min.regions = 0)
# h3k4me1_count <- as.matrix(cisTopicObject@count.matrix)
##use CreatGeneActivityMatrix.R to generate gene activity matrix


#RNA <- readRDS("../RNA_4cell.rds")
#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k4me1,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
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
write.table(h3k4me1_pre1, "../h3k4me1/h3k4me1_4cell_pseudo.txt", quote = F, sep = " ")

saveRDS(h3k4me1, "../h3k4me1/h3k4me1_4cell_integration.rds")
saveRDS(transfer.anchors, "../h3k4me1/h3k4me1_4cell_transfer_anchors.rds")

######################################### H3K27ac #######################################

setwd("../h3k27ac/")
set.seed(2)
h3k27ac_counts <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/h3k27ac.5kb.best.countMatrix.txt", sep =" ") 
h3k27ac_counts <- h3k27ac_counts[,grep("4cell",colnames(h3k27ac_counts))] 
h3k27ac <- CreateSeuratObject(counts = h3k27ac_counts,assay = 'peaks',project = 'h3k27ac',min.cells = 1)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = 'q0')
h3k27ac <- RunSVD(object = h3k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k27ac/h3k27ac_allStage_peaks_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(h3k27ac)]
h3k27ac[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k27ac) <- "ACTIVITY"
h3k27ac <- FindVariableFeatures(h3k27ac, nfeatures = 10000)
h3k27ac <- NormalizeData(h3k27ac)
h3k27ac <- ScaleData(h3k27ac)

#RNA <- readRDS("../RNA_4cell.rds")
#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k27ac,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
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
write.table(h3k27ac_pre1, "./h3k27ac_4cell_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k27ac, "h3k27ac_4cell_integration.rds")
saveRDS(transfer.anchors, "h3k27ac_4cell_transfer_anchors.rds")

############################################### H3K36me3 #########################################
#rm(list = ls())
setwd("../h3k36me3/")
set.seed(3)
load("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/h3k36me3.5kb.Rdata")
h3k36me3 <- subset(h3k36me3,cells = colnames(h3k36me3[,grep("4cell", h3k36me3@meta.data$group.ident)]))
pdf("h3k36me3_4cell_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 
####目前的UMAP不合适，有batch effect，会影响trajectory，所以重新降维
h3k36me3 <- RunUMAP(object = h3k36me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
h3k36me3 <- FindNeighbors(object = h3k36me3,reduction = 'lsi',dims = c(2:15))
h3k36me3 <- FindClusters(object = h3k36me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
pdf("h3k36me3_4cell_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

##直接去掉2241
h3k36me3 <- subset(h3k36me3,cells = colnames(h3k36me3[,grep("22101", h3k36me3@meta.data$tab.ident)]))
pdf("h3k36me3_4cell_only22101_5kbbin.pdf")
DepthCor(h3k36me3)
DimPlot(h3k36me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k36me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k36me3/h3k36me3_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(h3k36me3)]
h3k36me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k36me3) <- "ACTIVITY"
h3k36me3 <- FindVariableFeatures(h3k36me3, nfeatures = 10000)
h3k36me3 <- NormalizeData(h3k36me3)
h3k36me3 <- ScaleData(h3k36me3)

#RNA <- readRDS("../RNA_4cell.rds")
#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k36me3,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
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
write.table(h3k36me3_pre1, "./h3k36me3_4cell_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k36me3, "h3k36me3_4cell_integration.rds")
saveRDS(transfer.anchors, "h3k36me3_4cell_transfer_anchors.rds")

########################################### H3K27me3 #############################################
##借助coTACIT数据
setwd("../cotacit_data/")
pathToBams1 <- "./h3k27ac/"
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl5kb_all+.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='h3k27ac', min.regions = 0)
cotacit_count <- as.matrix(cisTopicObject@count.matrix)

cotacit <- CreateSeuratObject(counts = cotacit_count,assay = 'peaks',project = 'cotacit',min.cells = 1)
cotacit <- RunTFIDF(cotacit)
cotacit <- FindTopFeatures(cotacit, min.cutoff = 'q0')
cotacit <- RunSVD(object = cotacit, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

cotacit <- RunUMAP(object = cotacit,reduction = 'lsi',dims = c(2,4:25), n.neighbors = 10L)#或可调试reduction
cotacit <- FindNeighbors(object = cotacit,reduction = 'lsi',dims = c(2,4:25))
cotacit <- FindClusters(object = cotacit,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

label <- rownames(as.data.frame(cotacit@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("23011|23012",label[,1]),2] <- "morula"
label[grep("blasto|22115|22116",label[,1]),2] <- "blastocyst"

label[,3] <- cotacit@meta.data$seurat_clusters

label[grep("23021",label[,1]),4] <- "23021"
label[grep("23022",label[,1]),4] <- "23022"
label[grep("23023",label[,1]),4] <- "23023"
label[grep("23024",label[,1]),4] <- "23024"
label[grep("22115",label[,1]),4] <- "22115"
label[grep("22116",label[,1]),4] <- "22116"
label[grep("23011",label[,1]),4] <- "23011"
label[grep("23012",label[,1]),4] <- "23012"

cotacit@meta.data$group.ident <- as.factor(label[,2])
cotacit@meta.data$tab.ident <- as.factor(label[,4])

pdf("cotacit_all_5kbbin_h3k27ac.pdf")
DepthCor(cotacit)
DimPlot(cotacit, label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(cotacit, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(cotacit, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

saveRDS(cotacit, "cotacit_h3k27ac.rds")


cotacit_4cell <- subset(cotacit,cells = colnames(cotacit[,grep("4cell", cotacit@meta.data$group.ident)]))
geneActivity <- read.csv("./cotacit_h3k27ac_allstage_tss5kb_ActivityMatrix.txt", sep = " ",check.names=F)
geneActivity <- geneActivity[,colnames(cotacit_4cell)]
cotacit_4cell[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(cotacit_4cell) <- "ACTIVITY"
cotacit_4cell <- FindVariableFeatures(cotacit_4cell, nfeatures = 10000)
cotacit_4cell <- NormalizeData(cotacit_4cell)
cotacit_4cell <- ScaleData(cotacit_4cell)

#####integration#####
transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = cotacit_4cell,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = cotacit_4cell[["lsi"]], dims = 2:ncol(cotacit_4cell[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_pre1 <- ID.predict
length(unique(cotacit_pre1$predicted.id))
#cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]
write.table(cotacit_pre1, "./cotacit_4cell_h3k27ac_pseudo.txt", quote = F, sep = " ")
saveRDS(cotacit_4cell, "cotacit_4cell_h3k27ac_integration.rds")
saveRDS(transfer.anchors, "cotacit_4cell_h3k27ac_transfer_anchors.rds")


############################## H3K27me3 #################################################
setwd("../h3k27me3/")
set.seed(5)
#load("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/h3k27me3.Rdata")
h3k27me3 <- readRDS("./h3k27me3_15kbbin.rds")
h3k27me3 <- subset(h3k27me3,cells = colnames(h3k27me3[,grep("4cell", h3k27me3@meta.data$group.ident)]))
h3k27me3 <- RunUMAP(object = h3k27me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
h3k27me3 <- FindNeighbors(object = h3k27me3,reduction = 'lsi',dims = c(2:15))
h3k27me3 <- FindClusters(object = h3k27me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm
pdf("h3k27me3_4cell_5kbbin.pdf")
DepthCor(h3k27me3)
DimPlot(h3k27me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(h3k27me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

geneActivity <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/h3k27me3/h3k27me3_allstage_5kbbin_Activityrix.txt", sep = " ",check.names=F)
geneActivity <- 1/((geneActivity)+1)
geneActivity <- geneActivity[,colnames(h3k27me3)]
h3k27me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k27me3) <- "ACTIVITY"
h3k27me3 <- FindVariableFeatures(h3k27me3, nfeatures = 10000)
h3k27me3 <- NormalizeData(h3k27me3)
h3k27me3 <- ScaleData(h3k27me3)


transfer.anchors <- FindTransferAnchors(reference = RNA_pseudo,
                                        query = h3k27me3,
                                        features = VariableFeatures(object = RNA_pseudo),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5, 
                                        npcs = 10, dims = 1:10,
                                        k.filter = NA, 
                                        k.score = 10)
#Annotate scATAC-seq cells via label transfer
#transfer ID
RNA_pseudo@meta.data$ID <- colnames(RNA_pseudo)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA_pseudo$ID,
                                   weight.reduction = h3k27me3[["lsi"]], dims = 2:ncol(h3k27me3[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
h3k27me3_pre1 <- ID.predict
length(unique(h3k27me3_pre1$predicted.id))
#h3k27me3_pre1 <- h3k27me3_pre1[which(h3k27me3_pre1$prediction.score.max >0.2),]
write.table(h3k27me3_pre1, "./h3k27me3_4cell_pseudo.txt", quote = F, sep = " ")
saveRDS(h3k27me3, "h3k27me3_4cell_integration.rds")
saveRDS(transfer.anchors, "h3k27me3_4cell_transfer_anchors.rds")


####################################### all #######################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
h3k4me1_pre1 <- read.csv("./h3k4me1/h3k4me1_4cell_pseudo.txt", sep = " ")
h3k4me3_pre1 <- read.csv("./h3k4me3/h3k4me3_4cell_pseudo_inte.txt", sep = " ")
h3k27ac_pre1 <- read.csv("./h3k27ac/h3k27ac_4cell_pseudo.txt", sep = " ")
h3k36me3_pre1 <- read.csv("./h3k36me3/h3k36me3_4cell_pseudo.txt", sep = " ")
cotacit_pre1 <- read.csv("./cotacit_data/cotacit_4cell_h3k27ac_pseudo.txt", sep = " ")

h3k4me1_pre1 <- h3k4me1_pre1[which(h3k4me1_pre1$prediction.score.max >0.2),]
h3k4me3_pre1 <- h3k4me3_pre1[which(h3k4me3_pre1$prediction.score.max >0.2),]
h3k27ac_pre1 <- h3k27ac_pre1[which(h3k27ac_pre1$prediction.score.max >0.2),]
h3k36me3_pre1 <- h3k36me3_pre1[which(h3k36me3_pre1$prediction.score.max >0.2),]
cotacit_pre1 <- cotacit_pre1[which(cotacit_pre1$prediction.score.max >0.2),]

inte <- Reduce(intersect, list(h3k4me3_pre1$predicted.id, h3k4me1_pre1$predicted.id, 
                               h3k36me3_pre1$predicted.id, h3k27ac_pre1$predicted.id, cotacit_pre1$predicted.id))

h3k4me1_pre1 <- h3k4me1_pre1 [which(h3k4me1_pre1 $predicted.id %in% inte),]
h3k4me1_pre1 <- h3k4me1_pre1[order(h3k4me1_pre1$prediction.score.max, decreasing = T),]
h3k4me1_pre1  <- h3k4me1_pre1 [!duplicated(h3k4me1_pre1 $predicted.id),]
h3k4me1.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me1/allBam/Min/", rownames(h3k4me1_pre1), " ./4cell/", h3k4me1_pre1$predicted.id, "_inte.bam")
write.table(h3k4me1.sh, "./h3k4me1/h3k4me1_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1[order(h3k4me3_pre1$prediction.score.max, decreasing = T),]
h3k4me3_pre1  <- h3k4me3_pre1 [!duplicated(h3k4me3_pre1 $predicted.id),]
h3k4me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/", rownames(h3k4me3_pre1), " ./4cell/", h3k4me3_pre1$predicted.id, "_inte.bam")
write.table(h3k4me3.sh, "./h3k4me3/h3k4me3_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1[order(h3k27ac_pre1$prediction.score.max, decreasing = T),]
h3k27ac_pre1  <- h3k27ac_pre1 [!duplicated(h3k27ac_pre1 $predicted.id),]
h3k27ac.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/bam2000/", rownames(h3k27ac_pre1), " ./4cell/", h3k27ac_pre1$predicted.id, "_inte.bam")
write.table(h3k27ac.sh, "./h3k27ac/h3k27ac_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1[order(h3k36me3_pre1$prediction.score.max, decreasing = T),]
h3k36me3_pre1  <- h3k36me3_pre1 [!duplicated(h3k36me3_pre1 $predicted.id),]
h3k36me3.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/all_h3k36me3/bam2000/", rownames(h3k36me3_pre1), " ./4cell/", h3k36me3_pre1$predicted.id, "_inte.bam")
write.table(h3k36me3.sh, "./h3k36me3/h3k36me3_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

cotacit_pre1 <- cotacit_pre1 [which(cotacit_pre1 $predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1[order(cotacit_pre1$prediction.score.max, decreasing = T),]
cotacit_pre1  <- cotacit_pre1 [!duplicated(cotacit_pre1 $predicted.id),]
cotacit_pre1$ID <- rownames(cotacit_pre1)
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_all_h3k9me3_rmdup.bam")
cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2, " ./4cell/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k27me3/cotacit_4cell_linked_h3k27me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
cotacit.sh <- paste0("cp /media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3, " ./4cell/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k9me3/cotacit_4cell_linked_h3k9me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(cotacit_pre1[,c("predicted.id","prediction.score.max", "ID2", "ID3")], "./cotacit_data/integration/cotacit_4cell_integration.txt")



############################################# cotacit + tacit h3k37me3 #########################################

############ cotacit h3k27me3 ################
pathToBams1 <- "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/"
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/h3k27me3.peak-0.05.10kbMerge.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='h3k27ac', min.regions = 0)
cotacit_h3k27me3_count <- as.matrix(cisTopicObject@count.matrix)

cotacit_h3k27me3 <- CreateSeuratObject(counts = cotacit_h3k27me3_count,assay = 'peaks',project = 'cotacit_h3k27me3',min.cells = 1)
cotacit_h3k27me3 <- RunTFIDF(cotacit_h3k27me3)
cotacit_h3k27me3 <- FindTopFeatures(cotacit_h3k27me3, min.cutoff = 'q0')
cotacit_h3k27me3 <- RunSVD(object = cotacit_h3k27me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

cotacit_h3k27me3 <- RunUMAP(object = cotacit_h3k27me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
cotacit_h3k27me3 <- FindNeighbors(object = cotacit_h3k27me3,reduction = 'lsi',dims = c(2:15))
cotacit_h3k27me3 <- FindClusters(object = cotacit_h3k27me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

label <- rownames(as.data.frame(cotacit_h3k27me3@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("23011|23012",label[,1]),2] <- "morula"
label[grep("blasto|22115|22116",label[,1]),2] <- "blastocyst"

label[,3] <- cotacit_h3k27me3@meta.data$seurat_clusters

label[grep("23021",label[,1]),4] <- "23021"
label[grep("23022",label[,1]),4] <- "23022"
label[grep("23023",label[,1]),4] <- "23023"
label[grep("23024",label[,1]),4] <- "23024"
label[grep("22115",label[,1]),4] <- "22115"
label[grep("22116",label[,1]),4] <- "22116"
label[grep("23011",label[,1]),4] <- "23011"
label[grep("23012",label[,1]),4] <- "23012"

cotacit_h3k27me3@meta.data$group.ident <- as.factor(label[,2])
cotacit_h3k27me3@meta.data$tab.ident <- as.factor(label[,4])

pdf("./cotacit_data/cotacit_h3k27me3_peaks.pdf")
DepthCor(cotacit_h3k27me3)
DimPlot(cotacit_h3k27me3, label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(cotacit_h3k27me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(cotacit_h3k27me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

saveRDS(cotacit_h3k27me3, "./cotacit_data/cotacit_h3k27me3_peaks.rds")

#################### merge two datasets #######################
library(Signac)

pathToBams1 <- c('/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|ArchR|bai|bed",invert = T, bamFiles1)]
regions <- "/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/h3k27me3.peak-0.05.10kbMerge.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27me3')
h3k27me3_count <- as.matrix(cisTopicObject@count.matrix)

h3k27me3 <- CreateSeuratObject(counts = h3k27me3_count,assay = 'peaks',project = 'h3k27me3',min.cells = 1)
h3k27me3 <- RunTFIDF(h3k27me3)
h3k27me3 <- FindTopFeatures(h3k27me3, min.cutoff = 'q0')
h3k27me3 <- RunSVD(object = h3k27me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
h3k27me3 <- RunUMAP(object = h3k27me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
h3k27me3 <- FindNeighbors(object = h3k27me3,reduction = 'lsi',dims = c(2:15))
h3k27me3 <- FindClusters(object = h3k27me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

label <- rownames(as.data.frame(h3k27me3@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("E2.5",label[,1]),2] <- "morula"
label[grep("E3.5",label[,1]),2] <- "blastocyst"
h3k27me3@meta.data$group.ident <- as.factor(label[,2])

h3k27me3$dataset <- "tacit"
cotacit_h3k27me3$dataset <- "cotacit"

h3k27me3_4cell <- subset(h3k27me3,cells = colnames(h3k27me3)[grep("4cell", colnames(h3k27me3))])
cotacit_h3k27me3_4cell <- subset(cotacit_h3k27me3,cells = colnames(cotacit_h3k27me3)[grep("4cell", colnames(cotacit_h3k27me3))])
merge <- merge(h3k27me3_4cell, cotacit_h3k27me3_4cell)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_4cell.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k27me3_4cell, cotacit_h3k27me3_4cell),
  anchor.features = rownames(h3k27me3_4cell),
  reduction = "rlsi",
  dims = 2:30)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merge[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  k.weight=20)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k27me3_peaks_integration_4cell.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k27me3_4cell,
  query = cotacit_h3k27me3_4cell,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k27me3_4cell@meta.data$ID <- colnames(h3k27me3_4cell)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k27me3_4cell$ID,
                                   weight.reduction = cotacit_h3k27me3_4cell[["lsi"]], dims = 2:ncol(cotacit_h3k27me3_4cell[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k27me3_pre1 <- ID.predict
length(unique(cotacit_h3k27me3_pre1$predicted.id))

write.table(cotacit_h3k27me3_pre1, "./cotacit_data/h3k27me3_4cell_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k27me3_4cell_cotacit_tacit_inte_transfer_anchors.rds")
saveRDS(h3k27me3, "./cotacit_data/h3k27me3_peaks.rds")

############################################# cotacit + tacit h3k9me3 #########################################

############ cotacit h3k9me3 ################
pathToBams1 <- "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/"
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
regions <- "/media/helab/data1/min/00_reference/mm10.10K.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='h3k27ac', min.regions = 0)
cotacit_h3k9me3_count <- as.matrix(cisTopicObject@count.matrix)

cotacit_h3k9me3 <- CreateSeuratObject(counts = cotacit_h3k9me3_count,assay = 'peaks',project = 'cotacit_h3k9me3',min.cells = 1)
cotacit_h3k9me3 <- RunTFIDF(cotacit_h3k9me3)
cotacit_h3k9me3 <- FindTopFeatures(cotacit_h3k9me3, min.cutoff = 'q0')
cotacit_h3k9me3 <- RunSVD(object = cotacit_h3k9me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

cotacit_h3k9me3 <- RunUMAP(object = cotacit_h3k9me3,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
cotacit_h3k9me3 <- FindNeighbors(object = cotacit_h3k9me3,reduction = 'lsi',dims = c(2:15))
cotacit_h3k9me3 <- FindClusters(object = cotacit_h3k9me3,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

label <- rownames(as.data.frame(cotacit_h3k9me3@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("8cell",label[,1]),2] <- "8cell"
label[grep("23011|23012",label[,1]),2] <- "morula"
label[grep("blasto|22115|22116",label[,1]),2] <- "blastocyst"

label[,3] <- cotacit_h3k9me3@meta.data$seurat_clusters

label[grep("23021",label[,1]),4] <- "23021"
label[grep("23022",label[,1]),4] <- "23022"
label[grep("23023",label[,1]),4] <- "23023"
label[grep("23024",label[,1]),4] <- "23024"
label[grep("22115",label[,1]),4] <- "22115"
label[grep("22116",label[,1]),4] <- "22116"
label[grep("23011",label[,1]),4] <- "23011"
label[grep("23012",label[,1]),4] <- "23012"

cotacit_h3k9me3@meta.data$group.ident <- as.factor(label[,2])
cotacit_h3k9me3@meta.data$tab.ident <- as.factor(label[,4])

pdf("./cotacit_data/cotacit_h3k9me3_10kbbin.pdf")
DepthCor(cotacit_h3k9me3)
DimPlot(cotacit_h3k9me3, label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(cotacit_h3k9me3, group.by = "group.ident", label = TRUE, pt.size = 2,repel = TRUE)
DimPlot(cotacit_h3k9me3, group.by = "tab.ident", label = TRUE, pt.size = 2,repel = TRUE)
dev.off() 

saveRDS(cotacit_h3k9me3, "./cotacit_data/cotacit_h3k9me3_10kbbin.rds")

#################### merge two datasets #######################
library(Signac)
h3k9me3 <- readRDS("./cotacit_data/h3k9me3_10kbbin.rds")
h3k9me3$dataset <- "tacit"
cotacit_h3k9me3$dataset <- "cotacit"

h3k9me3_4cell <- subset(h3k9me3,cells = colnames(h3k9me3)[grep("4cell", colnames(h3k9me3))])
cotacit_h3k9me3_4cell <- subset(cotacit_h3k9me3,cells = colnames(cotacit_h3k9me3)[grep("4cell", colnames(cotacit_h3k9me3))])
merge <- merge(h3k9me3_4cell, cotacit_h3k9me3_4cell)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15), n.neighbors = 10L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.8,verbose = FALSE)#或可调试reduction和algorithm

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_4cell.pdf")
DimPlot(merge, reduction = "umap", label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "group.ident" ,label = TRUE, pt.size = 2)
DimPlot(object = merge, reduction = "umap",group.by = "dataset" ,label = TRUE, pt.size = 2)
dev.off()

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(h3k9me3_4cell, cotacit_h3k9me3_4cell),
  anchor.features = rownames(h3k9me3_4cell),
  reduction = "rlsi",
  dims = 2:30)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = merge[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  k.weight=20)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)

pdf("./cotacit_data/TACIT_coTACIT_merge_h3k9me3_peaks_integration_4cell.pdf")
DimPlot(integrated, group.by = "dataset")
DimPlot(integrated, group.by = "group.ident")
dev.off()

############ RNA imputation, transfer continuous data ################
# find transfer anchors
transfer.anchors <- FindTransferAnchors(
  reference = h3k9me3_4cell,
  query = cotacit_h3k9me3_4cell,
  reference.reduction = "lsi",
  reduction = "lsiproject",
  dims = 2:30
)

#transfer ID
h3k9me3_4cell@meta.data$ID <- colnames(h3k9me3_4cell)
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = h3k9me3_4cell$ID,
                                   weight.reduction = cotacit_h3k9me3_4cell[["lsi"]], dims = 2:ncol(cotacit_h3k9me3_4cell[["lsi"]]), k.weight = 10)
ID.predict <-  subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
cotacit_h3k9me3_pre1 <- ID.predict
length(unique(cotacit_h3k9me3_pre1$predicted.id))

write.table(cotacit_h3k9me3_pre1, "./cotacit_data/h3k9me3_4cell_cotacit_tacit_inte.txt", quote = F, sep = " ")
saveRDS(transfer.anchors, "./cotacit_data/h3k9me3_4cell_cotacit_tacit_inte_transfer_anchors.rds")

cotacit_pre1$k27_tacit <- cotacit_h3k27me3_pre1$predicted.id[match(cotacit_pre1$ID2, rownames(cotacit_h3k27me3_pre1))]
cotacit_pre1$k9_tacit <- cotacit_h3k9me3_pre1$predicted.id[match(cotacit_pre1$ID3, rownames(cotacit_h3k9me3_pre1))]
k27me3.sh <- paste0("samtools merge ./h3k27me3_tacit+cotacit/4cell/", cotacit_pre1$predicted.id, "_inte.bam ", 
                    "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k27me3/", cotacit_pre1$ID2,
                    " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k27me3/allBam/", cotacit_pre1$k27_tacit)
k9me3.sh <- paste0("samtools merge ./h3k9me3_tacit+cotacit/4cell/", cotacit_pre1$predicted.id, "_inte.bam ", 
                    "/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/cotacit_data/h3k9me3/", cotacit_pre1$ID3,
                    " /media/helab/data1/min/02_tacit/03_early_stage/all_h3k9me3/bam2000/", cotacit_pre1$k9_tacit)
write.table(k27me3.sh, "./cotacit_data/integration/h3k27me3_4cell_merge.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(k9me3.sh, "./cotacit_data/integration/h3k9me3_4cell_merge.sh", sep = "\t", quote = F, col.names = F, row.names = F)

write.table(cotacit_pre1, "./cotacit_data/cotacit_tacit_inte_4cell.txt", quote = F, sep = " ")



############################################### all ##############################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell")
merge.sh <- paste0("samtools merge ./all_agg_cor/", inte, "_all_marker.bam ", 
                   "../h3k4me1/4cell/", inte, "_inte.bam ", "../h3k4me3/4cell/", inte, "_inte.bam ",
                   "../h3k27ac/4cell/", inte, "_inte.bam ", "../h3k36me3/4cell/", inte, "_inte.bam ",
                   "../cotacit_data/integration/h3k9me3_tacit+cotacit/4cell/", inte, "_inte.bam ",
                   "../cotacit_data/integration/h3k27me3_tacit+cotacit/4cell/", inte, "_inte.bam ")
write.table(merge.sh, "merge_all_markers.sh", sep = "\t", quote = F, col.names = F, row.names = F)


############################### plote for chromHMM ##############################################
library(pheatmap)
C2 <- read.csv("./4cell_4/learnmodel/emissions_12_chrhmm.txt", sep="\t")
rownames(C2) <- C2[,1]
C2 <- C2[,-1]
dat <- as.matrix(C2)
colnames(dat) <- colnames(C2)
rownames(dat) <- rownames(C2)
dat[dat<=0.01] <- 0
dat <- log10(dat+1)
pdf("4cell_4_intersect_nochrM_emissions_10_chrhmm.pdf")
p <- pheatmap(dat, 
              cluster_rows =F, cluster_cols = F, 
              show_rownames = T, scale ="none",
              show_colnames = T,
              clustering_method="ward.D2",
              clustering_distance_rows = "correlation",
              color = colorRampPalette(c("#D2EAF9","#277FC1","#2172B7"))(500))
#clustering_method should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
p
dev.off()

############################## annotation chromatin states ########################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
######## 4cell_1 #########
C2_bed <- read.csv("./03_chromHMM_4cell/4cell_1/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "K27"
C2_bed[which(C2_bed$V4 == 2),10] <- "He"
C2_bed[which(C2_bed$V4 == 3),10] <- "Un"
C2_bed[which(C2_bed$V4 == 4),10] <- "K9"
C2_bed[which(C2_bed$V4 == 5),10] <- "He"
C2_bed[which(C2_bed$V4 == 6),10] <- "K9"
C2_bed[which(C2_bed$V4 == 7),10] <- "Un"
C2_bed[which(C2_bed$V4 == 8),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 9),10] <- "Es"
C2_bed[which(C2_bed$V4 == 10),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 11),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ew"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "03_chromHMM_4cell/4cell_1/4cell_1_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")


C2_bed_E <- C2_bed[grep("Ep|Es|Ew", C2_bed$V10),]
write.table(C2_bed_E, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_enhancer.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_P <- C2_bed[grep("Pw|Ps", C2_bed$V10),]
write.table(C2_bed_P, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_promoter.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_T <- C2_bed[grep("Tw|Ts", C2_bed$V10),]
write.table(C2_bed_T, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_transcription.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_M <- C2_bed[grep("M", C2_bed$V10),]
write.table(C2_bed_M, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_multivalent.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K9 <- C2_bed[grep("K9", C2_bed$V10),]
write.table(C2_bed_K9, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_h3k9me3.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("K27", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_polycomb.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("He", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_1/4cell_1_chrhmm_heterochromtin.bed", quote = F, col.names = F, row.names = F, sep = "\t")

####### change colors for bed file ##############
C1_bed <- read.csv("03_chromHMM_4cell/4cell_1/4cell_1_mm10_chrhmm_dense_anno.bed", sep = "\t", header = F)
C1_bed[grep("Un", C1_bed$V4),9] <- "209,211,212"
C1_bed[grep("Mu", C1_bed$V4),9] <- "170,104,55"
C1_bed[grep("Pw", C1_bed$V4),9] <- "117,204,158"
C1_bed[grep("Ps", C1_bed$V4),9] <- "11,158,91"
C1_bed[grep("Ew", C1_bed$V4),9] <- "249,209,160"
C1_bed[grep("Es", C1_bed$V4),9] <- "249,156,8"
C1_bed[grep("Tw", C1_bed$V4),9] <- "214,107,123"
C1_bed[grep("Ts", C1_bed$V4),9] <- "198,32,48"
C1_bed[grep("K27", C1_bed$V4),9] <- "146,234,244"
C1_bed[grep("K9", C1_bed$V4),9] <- "21,162,206"
C1_bed[grep("He", C1_bed$V4),9] <- "109,110,113"
write.table(C1_bed, "03_chromHMM_4cell/4cell_1/4cell_1_mm10_chrhmm_dense_anno_color.bed", quote = F, col.names = F, row.names = F, sep = "\t")


######### 4cell_2 ##########
C2_bed <- read.csv("./03_chromHMM_4cell/4cell_2/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Es"
C2_bed[which(C2_bed$V4 == 2),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 3),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 4),10] <- "Ew"
C2_bed[which(C2_bed$V4 == 5),10] <- "Un"
C2_bed[which(C2_bed$V4 == 6),10] <- "Un"
C2_bed[which(C2_bed$V4 == 7),10] <- "He"
C2_bed[which(C2_bed$V4 == 8),10] <- "K9"
C2_bed[which(C2_bed$V4 == 9),10] <- "He"
C2_bed[which(C2_bed$V4 == 10),10] <- "K27"
C2_bed[which(C2_bed$V4 == 11),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 12),10] <- "Ps"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "03_chromHMM_4cell/4cell_2/4cell_2_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")


C2_bed_E <- C2_bed[grep("Ep|Es|Ew", C2_bed$V10),]
write.table(C2_bed_E, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_enhancer.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_P <- C2_bed[grep("Pw|Ps", C2_bed$V10),]
write.table(C2_bed_P, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_promoter.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_T <- C2_bed[grep("Tw|Ts", C2_bed$V10),]
write.table(C2_bed_T, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_transcription.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_M <- C2_bed[grep("M", C2_bed$V10),]
write.table(C2_bed_M, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_multivalent.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K9 <- C2_bed[grep("K9", C2_bed$V10),]
write.table(C2_bed_K9, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_h3k9me3.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("K27", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_polycomb.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("He", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_2/4cell_2_chrhmm_heterochromtin.bed", quote = F, col.names = F, row.names = F, sep = "\t")

####### change colors for bed file ##############
C1_bed <- read.csv("03_chromHMM_4cell/4cell_2/4cell_2_mm10_chrhmm_dense_anno.bed", sep = "\t", header = F)
C1_bed[grep("Un", C1_bed$V4),9] <- "209,211,212"
C1_bed[grep("Mu", C1_bed$V4),9] <- "170,104,55"
C1_bed[grep("Pw", C1_bed$V4),9] <- "117,204,158"
C1_bed[grep("Ps", C1_bed$V4),9] <- "11,158,91"
C1_bed[grep("Ew", C1_bed$V4),9] <- "249,209,160"
C1_bed[grep("Es", C1_bed$V4),9] <- "249,156,8"
C1_bed[grep("Tw", C1_bed$V4),9] <- "214,107,123"
C1_bed[grep("Ts", C1_bed$V4),9] <- "198,32,48"
C1_bed[grep("K27", C1_bed$V4),9] <- "146,234,244"
C1_bed[grep("K9", C1_bed$V4),9] <- "21,162,206"
C1_bed[grep("He", C1_bed$V4),9] <- "109,110,113"
write.table(C1_bed, "03_chromHMM_4cell/4cell_2/4cell_2_mm10_chrhmm_dense_anno_color.bed", quote = F, col.names = F, row.names = F, sep = "\t")


######### 4cell_3 ##########
C2_bed <- read.csv("./03_chromHMM_4cell/4cell_3/learnmodel_inte/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Pw"
C2_bed[which(C2_bed$V4 == 2),10] <- "He"
C2_bed[which(C2_bed$V4 == 3),10] <- "He"
C2_bed[which(C2_bed$V4 == 4),10] <- "K27"
C2_bed[which(C2_bed$V4 == 5),10] <- "Un"
C2_bed[which(C2_bed$V4 == 6),10] <- "K9"
C2_bed[which(C2_bed$V4 == 7),10] <- "Un"
C2_bed[which(C2_bed$V4 == 8),10] <- "Un"
C2_bed[which(C2_bed$V4 == 9),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 10),10] <- "Tw"
C2_bed[which(C2_bed$V4 == 11),10] <- "Es"
C2_bed[which(C2_bed$V4 == 12),10] <- "Es"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "03_chromHMM_4cell/4cell_3/4cell_3_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")


C2_bed_E <- C2_bed[grep("Ep|Es|Ew", C2_bed$V10),]
write.table(C2_bed_E, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_enhancer.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_P <- C2_bed[grep("Pw|Ps", C2_bed$V10),]
write.table(C2_bed_P, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_promoter.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_T <- C2_bed[grep("Tw|Ts", C2_bed$V10),]
write.table(C2_bed_T, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_transcription.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_M <- C2_bed[grep("M", C2_bed$V10),]
write.table(C2_bed_M, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_multivalent.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K9 <- C2_bed[grep("K9", C2_bed$V10),]
write.table(C2_bed_K9, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_h3k9me3.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("K27", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_polycomb.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("He", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_3/4cell_3_chrhmm_heterochromtin.bed", quote = F, col.names = F, row.names = F, sep = "\t")

####### change colors for bed file ##############
C1_bed <- read.csv("03_chromHMM_4cell/4cell_3/4cell_3_mm10_chrhmm_dense_anno.bed", sep = "\t", header = F)
C1_bed[grep("Un", C1_bed$V4),9] <- "209,211,212"
C1_bed[grep("Mu", C1_bed$V4),9] <- "170,104,55"
C1_bed[grep("Pw", C1_bed$V4),9] <- "117,204,158"
C1_bed[grep("Ps", C1_bed$V4),9] <- "11,158,91"
C1_bed[grep("Ew", C1_bed$V4),9] <- "249,209,160"
C1_bed[grep("Es", C1_bed$V4),9] <- "249,156,8"
C1_bed[grep("Tw", C1_bed$V4),9] <- "214,107,123"
C1_bed[grep("Ts", C1_bed$V4),9] <- "198,32,48"
C1_bed[grep("K27", C1_bed$V4),9] <- "146,234,244"
C1_bed[grep("K9", C1_bed$V4),9] <- "21,162,206"
C1_bed[grep("He", C1_bed$V4),9] <- "109,110,113"
write.table(C1_bed, "03_chromHMM_4cell/4cell_3/4cell_3_mm10_chrhmm_dense_anno_color.bed", quote = F, col.names = F, row.names = F, sep = "\t")


######### 4cell_4 ##########
C2_bed <- read.csv("./03_chromHMM_4cell/4cell_4/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]
C2_bed[which(C2_bed$V4 == 1),10] <- "Es"
C2_bed[which(C2_bed$V4 == 2),10] <- "Es"
C2_bed[which(C2_bed$V4 == 3),10] <- "Ps"
C2_bed[which(C2_bed$V4 == 4),10] <- "Es"
C2_bed[which(C2_bed$V4 == 5),10] <- "Pw"
C2_bed[which(C2_bed$V4 == 6),10] <- "He"
C2_bed[which(C2_bed$V4 == 7),10] <- "He"
C2_bed[which(C2_bed$V4 == 8),10] <- "Un"
C2_bed[which(C2_bed$V4 == 9),10] <- "Un"
C2_bed[which(C2_bed$V4 == 10),10] <- "Ts"
C2_bed[which(C2_bed$V4 == 11),10] <- "Un"
C2_bed[which(C2_bed$V4 == 12),10] <- "He"
C2_bed <- C2_bed[,c("V1", "V2", "V3", "V10","V5", "V6", "V7", "V8", "V9")]
write.table(C2_bed, "03_chromHMM_4cell/4cell_4/4cell_4_mm10_chrhmm_dense_anno.bed", quote = F, col.names = F, row.names = F, sep = "\t")


C2_bed_E <- C2_bed[grep("Ep|Es|Ew", C2_bed$V10),]
write.table(C2_bed_E, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_enhancer.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_P <- C2_bed[grep("Pw|Ps", C2_bed$V10),]
write.table(C2_bed_P, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_promoter.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_T <- C2_bed[grep("Tw|Ts", C2_bed$V10),]
write.table(C2_bed_T, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_transcription.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_M <- C2_bed[grep("M", C2_bed$V10),]
write.table(C2_bed_M, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_multivalent.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K9 <- C2_bed[grep("K9", C2_bed$V10),]
write.table(C2_bed_K9, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_h3k9me3.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("K27", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_polycomb.bed", quote = F, col.names = F, row.names = F, sep = "\t")
C2_bed_K27 <- C2_bed[grep("He", C2_bed$V10),]
write.table(C2_bed_K27, "03_chromHMM_4cell/4cell_4/4cell_4_chrhmm_heterochromtin.bed", quote = F, col.names = F, row.names = F, sep = "\t")

####### change colors for bed file ##############
C1_bed <- read.csv("03_chromHMM_4cell/4cell_4/4cell_4_mm10_chrhmm_dense_anno.bed", sep = "\t", header = F)
C1_bed[grep("Un", C1_bed$V4),9] <- "209,211,212"
C1_bed[grep("Mu", C1_bed$V4),9] <- "170,104,55"
C1_bed[grep("Pw", C1_bed$V4),9] <- "117,204,158"
C1_bed[grep("Ps", C1_bed$V4),9] <- "11,158,91"
C1_bed[grep("Ew", C1_bed$V4),9] <- "249,209,160"
C1_bed[grep("Es", C1_bed$V4),9] <- "249,156,8"
C1_bed[grep("Tw", C1_bed$V4),9] <- "214,107,123"
C1_bed[grep("Ts", C1_bed$V4),9] <- "198,32,48"
C1_bed[grep("K27", C1_bed$V4),9] <- "146,234,244"
C1_bed[grep("K9", C1_bed$V4),9] <- "21,162,206"
C1_bed[grep("He", C1_bed$V4),9] <- "109,110,113"
write.table(C1_bed, "03_chromHMM_4cell/4cell_4/4cell_4_mm10_chrhmm_dense_anno_color.bed", quote = F, col.names = F, row.names = F, sep = "\t")


################################# zga and toti gene expression ############################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration")
RNA_pseudo <- readRDS("scRNA/RNA_pseudo.rds")
RNA_pseudo <- as.data.frame(RNA_pseudo@assays[["RNA"]]@counts)
RNA_pseudo <- RNA_pseudo[,grep("zygote|2cell|4cell", colnames(RNA_pseudo))]

toti <- read.csv("./totipotency.markers.csv", sep = ",", header = F)
pluri<- read.csv("./pluripotency.markers.csv", sep = ",", header = F)
ZGA <- read.xlsx("/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/MajorZGA_genes.xlsx")
ZGA <- na.omit(ZGA$Gene)
RNA_ZGA <- RNA_pseudo[which(rownames(RNA_pseudo) %in% ZGA),]
RNA_toti <- RNA_pseudo[which(rownames(RNA_pseudo) %in% toti$V1),]
RNA_pluri <- RNA_pseudo[which(rownames(RNA_pseudo) %in% pluri$V1),]
RNA_ZGA <- colSums(RNA_ZGA)/620
RNA_toti <- colSums(RNA_toti)/215
RNA_pluri <- colSums(RNA_pluri)/38

dat <- data.frame(cells=colnames(RNA_pseudo), zga=RNA_ZGA, toti=RNA_toti, pluri=RNA_pluri)
dat[grep("zygote", dat$cells),5] <- "zygote"
dat[which(dat$cells == "2cell_pseudo_1" | dat$cells == "2cell_pseudo_6" | dat$cells == "2cell_pseudo_8" |dat$cells == "2cell_pseudo_9" |dat$cells == "2cell_pseudo_13"),5] <- "2cell_1"
dat[which(dat$cells == "2cell_pseudo_14" | dat$cells == "2cell_pseudo_22" | dat$cells == "2cell_pseudo_24" |dat$cells == "2cell_pseudo_25" |dat$cells == "2cell_pseudo_26"),5] <- "2cell_2"
dat[which(dat$cells == "4cell_pseudo_1" | dat$cells == "4cell_pseudo_2" | dat$cells == "4cell_pseudo_3" |dat$cells == "4cell_pseudo_5" |dat$cells == "4cell_pseudo_6"),5] <- "4cell_1"
dat[which(dat$cells == "4cell_pseudo_9" | dat$cells == "4cell_pseudo_11" | dat$cells == "4cell_pseudo_13" |dat$cells == "4cell_pseudo_14" ),5] <- "4cell_2"
dat[which(dat$cells == "4cell_pseudo_15" | dat$cells == "4cell_pseudo_16" | dat$cells == "4cell_pseudo_17" |dat$cells == "4cell_pseudo_18" ),5] <- "4cell_3"
dat[which(dat$cells == "4cell_pseudo_19" | dat$cells == "4cell_pseudo_20" | dat$cells == "4cell_pseudo_22" |dat$cells == "4cell_pseudo_23"),5] <- "4cell_4"
dat <- na.omit(dat)

pdf("pseudo_zygote_2cell_4cell_zga_toti.pdf")
ggplot(data=dat,mapping = aes(x = V5, y = pluri))+
  #geom_violin(aes(fill = V5), trim = FALSE) + 
  geom_boxplot(width = 0.5)+
  #scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0"))+
  #theme(legend.position = "none")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
dev.off()

############## use single cells #############
RNA <- readRDS("./scRNA/all_RNA_no_outlier.rds")
RNA_count <- as.data.frame(RNA@assays[["RNA"]]@counts)
RNA_count <- RNA_count[,grep("zygote|2cell|4cell", colnames(RNA_count))]

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

dat <- data.frame(cells=colnames(RNA_count), zga=RNA_ZGA, toti=RNA_toti, pluri=RNA_pluri)
mat <- read.csv("./scRNA/RNA_pseudo_id.txt", sep = "\t")

dat$pseudo_id <- mat$pseudo_id[match(dat$cells, mat$X)]
dat[grep("zygote", dat$pseudo_id),6] <- "zygote"
dat[which(dat$pseudo_id == "2cell_pseudo_1" | dat$pseudo_id == "2cell_pseudo_6" | dat$pseudo_id == "2cell_pseudo_8" |dat$pseudo_id == "2cell_pseudo_9" |dat$pseudo_id == "2cell_pseudo_13"),6] <- "2cell_1"
dat[which(dat$pseudo_id == "2cell_pseudo_14" | dat$pseudo_id == "2cell_pseudo_22" | dat$pseudo_id == "2cell_pseudo_24" |dat$pseudo_id == "2cell_pseudo_25" |dat$pseudo_id == "2cell_pseudo_26"),6] <- "2cell_2"
dat[which(dat$pseudo_id == "4cell_pseudo_1" | dat$pseudo_id == "4cell_pseudo_2" | dat$pseudo_id == "4cell_pseudo_3" |dat$pseudo_id == "4cell_pseudo_5" |dat$pseudo_id == "4cell_pseudo_6"),6] <- "4cell_1"
dat[which(dat$pseudo_id == "4cell_pseudo_9" | dat$pseudo_id == "4cell_pseudo_11" | dat$pseudo_id == "4cell_pseudo_13" |dat$pseudo_id == "4cell_pseudo_14" ),6] <- "4cell_2"
dat[which(dat$pseudo_id == "4cell_pseudo_15" | dat$pseudo_id == "4cell_pseudo_16" | dat$pseudo_id == "4cell_pseudo_17" |dat$pseudo_id == "4cell_pseudo_18" ),6] <- "4cell_3"
dat[which(dat$pseudo_id == "4cell_pseudo_19" | dat$pseudo_id == "4cell_pseudo_20" | dat$pseudo_id == "4cell_pseudo_22" |dat$pseudo_id == "4cell_pseudo_23"),6] <- "4cell_4"
dat <- na.omit(dat)

pdf("zygote_2cell_4cell_zga_toti.pdf")
ggplot(data=dat,mapping = aes(x = V6, y = pluri))+
  #geom_violin(aes(fill = V6), trim = FALSE) + 
  geom_boxplot(width = 0.5)+
  #scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0"))+
  #theme(legend.position = "none")+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
dev.off()

################################## state expression ###########################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell")

C2_bed <- read.csv("./4cell_4/learnmodel/mm10_12_chrhmm_dense.bed", sep = "\t", header = F)
C2_bed <- C2_bed[-1,]

E1 <- C2_bed[which(C2_bed$V4 == 1),]
E2 <- C2_bed[which(C2_bed$V4 == 2),]
E3 <- C2_bed[which(C2_bed$V4 == 3),]
E4 <- C2_bed[which(C2_bed$V4 == 4),]
E5 <- C2_bed[which(C2_bed$V4 == 5),]
E6 <- C2_bed[which(C2_bed$V4 == 6),]
E7 <- C2_bed[which(C2_bed$V4 == 7),]
E8 <- C2_bed[which(C2_bed$V4 == 8),]
E9 <- C2_bed[which(C2_bed$V4 == 9),]
E10 <- C2_bed[which(C2_bed$V4 == 10),]
E11 <- C2_bed[which(C2_bed$V4 == 11),]
E12 <- C2_bed[which(C2_bed$V4 == 12),]
write.table(E1, "./4cell_4/4cell_4_mm10_chrhmm_E1.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E2, "./4cell_4/4cell_4_mm10_chrhmm_E2.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E3, "./4cell_4/4cell_4_mm10_chrhmm_E3.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E4, "./4cell_4/4cell_4_mm10_chrhmm_E4.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E5, "./4cell_4/4cell_4_mm10_chrhmm_E5.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E6, "./4cell_4/4cell_4_mm10_chrhmm_E6.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E7, "./4cell_4/4cell_4_mm10_chrhmm_E7.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E8, "./4cell_4/4cell_4_mm10_chrhmm_E8.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E9, "./4cell_4/4cell_4_mm10_chrhmm_E9.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E10, "./4cell_4/4cell_4_mm10_chrhmm_E10.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E11, "./4cell_4/4cell_4_mm10_chrhmm_E11.bed", quote = F, col.names = F, row.names = F, sep = "\t")
write.table(E12, "./4cell_4/4cell_4_mm10_chrhmm_E12.bed", quote = F, col.names = F, row.names = F, sep = "\t")

###### annotate bed with Chipseeker #####
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

path <- "./4cell_4/"
files <- paste(path, list.files(path), sep='')
files <- files[grep(".bed", files)]
for (i in 1:12) {
  peak.anno <- annotatePeak(files[[i]], tssRegion=c(-3000, 3000),
                            TxDb=txdb, annoDb="org.Mm.eg.db")
  anno <- as.data.frame(peak.anno)
  write.table(anno, paste0(files[[i]], "_anno.txt"), quote = F, col.names = T, row.names = F, sep = "\t")
}


##################### expression for each chromatin state ###############################
RNA <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/all_RNA_no_outlier.rds")
RNA_count <- as.data.frame(RNA@assays[["RNA"]]@counts)
RNA_count <- RNA_count[,grep("zygote|2cell|4cell|8cell", colnames(RNA_count))]

dat <- data.frame(cells=colnames(RNA_count))
mat <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/RNA_pseudo_id.txt", sep = "\t")
dat$pseudo_id <- mat$pseudo_id[match(dat$cells, mat$X)]
dat[grep("zygote", dat$pseudo_id),3] <- "zygote"
dat[grep("8cell", dat$pseudo_id),3] <- "8cell"
dat[which(dat$pseudo_id == "2cell_pseudo_1" | dat$pseudo_id == "2cell_pseudo_6" | dat$pseudo_id == "2cell_pseudo_8" |dat$pseudo_id == "2cell_pseudo_9" |dat$pseudo_id == "2cell_pseudo_13"),3] <- "2cell_1"
dat[which(dat$pseudo_id == "2cell_pseudo_14" | dat$pseudo_id == "2cell_pseudo_22" | dat$pseudo_id == "2cell_pseudo_24" |dat$pseudo_id == "2cell_pseudo_25" |dat$pseudo_id == "2cell_pseudo_26"),3] <- "2cell_2"
dat[which(dat$pseudo_id == "4cell_pseudo_1" | dat$pseudo_id == "4cell_pseudo_2" | dat$pseudo_id == "4cell_pseudo_3" |dat$pseudo_id == "4cell_pseudo_5" |dat$pseudo_id == "4cell_pseudo_6"),3] <- "4cell_1"
dat[which(dat$pseudo_id == "4cell_pseudo_9" | dat$pseudo_id == "4cell_pseudo_11" | dat$pseudo_id == "4cell_pseudo_13" |dat$pseudo_id == "4cell_pseudo_14" ),3] <- "4cell_2"
dat[which(dat$pseudo_id == "4cell_pseudo_15" | dat$pseudo_id == "4cell_pseudo_16" | dat$pseudo_id == "4cell_pseudo_17" |dat$pseudo_id == "4cell_pseudo_18" ),3] <- "4cell_3"
dat[which(dat$pseudo_id == "4cell_pseudo_19" | dat$pseudo_id == "4cell_pseudo_20" | dat$pseudo_id == "4cell_pseudo_22" |dat$pseudo_id == "4cell_pseudo_23"),3] <- "4cell_4"

id_2cell_c1 <- dat[grep("2cell_1", dat$V3),1]
id_2cell_c2 <- dat[grep("2cell_2", dat$V3),1]
id_4cell_c1 <- dat[grep("4cell_1", dat$V3),1]
id_4cell_c2 <- dat[grep("4cell_2", dat$V3),1]
id_4cell_c3 <- dat[grep("4cell_3", dat$V3),1]
id_4cell_c4 <- dat[grep("4cell_4", dat$V3),1]

RNA_count_zygote <- RNA_count[,grep("zygote", colnames(RNA_count))]
RNA_count_8cell <- RNA_count[,grep("8cell", colnames(RNA_count))]
RNA_count_2cellC1 <- RNA_count[,which(colnames(RNA_count) %in% id_2cell_c1)]
RNA_count_2cellC2 <- RNA_count[,which(colnames(RNA_count) %in% id_2cell_c2)]
RNA_count_4cellC1 <- RNA_count[,which(colnames(RNA_count) %in% id_4cell_c1)]
RNA_count_4cellC2 <- RNA_count[,which(colnames(RNA_count) %in% id_4cell_c2)]
RNA_count_4cellC3 <- RNA_count[,which(colnames(RNA_count) %in% id_4cell_c3)]
RNA_count_4cellC4 <- RNA_count[,which(colnames(RNA_count) %in% id_4cell_c4)]

exp <- data.frame(zygote=rowMeans(RNA_count_zygote), C2_1=rowMeans(RNA_count_2cellC1), C2_2=rowMeans(RNA_count_2cellC2),
                  C4_1=rowMeans(RNA_count_4cellC1), C4_2=rowMeans(RNA_count_4cellC2), C4_3=rowMeans(RNA_count_4cellC3), C4_4=rowMeans(RNA_count_4cellC4),
                  C8_all=rowMeans(RNA_count_8cell))


##### 4cell_4 ###
E1 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E1.bed_anno.txt", sep = "\t", header = T)
E1 <- E1$SYMBOL
E1 <- E1[2:length(E1)]
E1 <- unique(E1)
c4_1_E1 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E1),]
c4_1_E1 <- rowMeans(c4_1_E1)
c4_1_E1 <- log2(c4_1_E1+1)
c4_1_E1 <- data.frame(state="E1", exp=c4_1_E1)


E2 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E2.bed_anno.txt", sep = "\t", header = T)
E2 <- E2$SYMBOL
E2 <- E2[2:length(E2)]
E2 <- unique(E2)
c4_1_E2 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E2),]
c4_1_E2 <- rowMeans(c4_1_E2)
c4_1_E2 <- log2(c4_1_E2+1)
c4_1_E2 <- data.frame(state="E2", exp=c4_1_E2)

E3 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E3.bed_anno.txt", sep = "\t", header = T)
E3 <- E3$SYMBOL
E3 <- E3[2:length(E3)]
E3 <- unique(E3)
c4_1_E3 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E3),]
c4_1_E3 <- rowMeans(c4_1_E3)
c4_1_E3 <- log2(c4_1_E3+1)
c4_1_E3 <- data.frame(state="E3", exp=c4_1_E3)

E4 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E4.bed_anno.txt", sep = "\t", header = T)
E4 <- E4$SYMBOL
E4 <- E4[2:length(E4)]
E4 <- unique(E4)
c4_1_E4 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E4),]
c4_1_E4 <- rowMeans(c4_1_E4)
c4_1_E4 <- log2(c4_1_E4+1)
c4_1_E4 <- data.frame(state="E4", exp=c4_1_E4)

E5 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E5.bed_anno.txt", sep = "\t", header = T)
E5 <- E5$SYMBOL
E5 <- E5[2:length(E5)]
E5 <- unique(E5)
c4_1_E5 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E5),]
c4_1_E5 <- rowMeans(c4_1_E5)
c4_1_E5 <- log2(c4_1_E5+1)
c4_1_E5 <- data.frame(state="E5", exp=c4_1_E5)

E6 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E6.bed_anno.txt", sep = "\t", header = T)
E6 <- E6$SYMBOL
E6 <- E6[2:length(E6)]
E6 <- unique(E6)
c4_1_E6 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E6),]
c4_1_E6 <- rowMeans(c4_1_E6)
c4_1_E6 <- log2(c4_1_E6+1)
c4_1_E6 <- data.frame(state="E6", exp=c4_1_E6)

E7 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E7.bed_anno.txt", sep = "\t", header = T)
E7 <- E7$SYMBOL
E7 <- E7[2:length(E7)]
E7 <- unique(E7)
c4_1_E7 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E7),]
c4_1_E7 <- rowMeans(c4_1_E7)
c4_1_E7 <- log2(c4_1_E7+1)
c4_1_E7 <- data.frame(state="E7", exp=c4_1_E7)

E8 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E8.bed_anno.txt", sep = "\t", header = T)
E8 <- E8$SYMBOL
E8 <- E8[2:length(E8)]
E8 <- unique(E8)
c4_1_E8 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E8),]
c4_1_E8 <- rowMeans(c4_1_E8)
c4_1_E8 <- log2(c4_1_E8+1)
c4_1_E8 <- data.frame(state="E8", exp=c4_1_E8)

E9 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E9.bed_anno.txt", sep = "\t", header = T)
E9 <- E9$SYMBOL
E9 <- E9[2:length(E9)]
E9 <- unique(E9)
c4_1_E9 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E9),]
c4_1_E9 <- rowMeans(c4_1_E9)
c4_1_E9 <- log2(c4_1_E9+1)
c4_1_E9 <- data.frame(state="E9", exp=c4_1_E9)

E10 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E10.bed_anno.txt", sep = "\t", header = T)
E10 <- E10$SYMBOL
E10 <- E10[2:length(E10)]
E10 <- unique(E10)
c4_1_E10 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E10),]
c4_1_E10 <- rowMeans(c4_1_E10)
c4_1_E10 <- log2(c4_1_E10+1)
c4_1_E10 <- data.frame(state="E10", exp=c4_1_E10)

E11 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E11.bed_anno.txt", sep = "\t", header = T)
E11 <- E11$SYMBOL
E11 <- E11[2:length(E11)]
E11 <- unique(E11)
c4_1_E11 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E11),]
c4_1_E11 <- rowMeans(c4_1_E11)
c4_1_E11 <- log2(c4_1_E11+1)
c4_1_E11 <- data.frame(state="E11", exp=c4_1_E11)

E12 <- read.csv("./4cell_4/4cell_4_mm10_chrhmm_E12.bed_anno.txt", sep = "\t", header = T)
E12 <- E12$SYMBOL
E12 <- E12[2:length(E12)]
E12 <- unique(E12)
c4_1_E12 <- RNA_count_4cellC4[which(rownames(RNA_count_4cellC4) %in% E12),]
c4_1_E12 <- rowMeans(c4_1_E12)
c4_1_E12 <- log2(c4_1_E12+1)
c4_1_E12 <- data.frame(state="E12", exp=c4_1_E12)

state_exp <- rbind(c4_1_E1, c4_1_E2, c4_1_E3, c4_1_E4, c4_1_E5, c4_1_E6, c4_1_E7, c4_1_E8, c4_1_E9, c4_1_E10, c4_1_E11, c4_1_E12)
pdf("./4cell_4/4cell_4_all_state_exp.pdf")
ggplot(data=state_exp,mapping = aes(x = state, y = exp))+
  #geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(aes(fill = state),width = 0.6)+
  #ylim(0,14)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=12),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())  #不显示网格线
dev.off()

write.table(state_exp, "./4cell_4/4cell_4_all_state_expression.txt", quote = F, sep = "\t")


################################### overlap with genome features (supplementary figure 4c)############################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell")
c4_1 <- read.csv("./4cell_1/learnmodel/mm10_12_chrhmm_overlap.txt", sep="\t")
c4_1 <- c4_1[,2:8]
c4_1 <- c4_1[c(12,9,8,10,11,5,2,1,6,4,3,7),]##与figure2c中的states顺序对应好，图可直接用
c4_1 <- c4_1[,-1]
library(pheatmap)
#c2_1 <- apply(c2_1, 2, function(x) log10(x+0.1))
pdf("./c4_1_state_overlap.pdf", width = 4, height = 6)
pheatmap(c4_1, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = T, scale ="row",
         show_colnames = T,
         color = colorRampPalette(c("#FFFFFF","#5050f9"))(100))
dev.off()


c4_2 <- read.csv("./4cell_2/learnmodel/mm10_12_chrhmm_overlap.txt", sep="\t")
c4_2 <- c4_2[,2:8]
c4_2 <- c4_2[c(2,1,3,4,11,12,9,7,10,8,5,6),]##与figure2c中的states顺序对应好，图可直接用
c4_2 <- c4_2[,-1]
library(pheatmap)
#c2_1 <- apply(c2_1, 2, function(x) log10(x+0.1))
pdf("./c4_2_state_overlap.pdf", width = 4, height = 6)
pheatmap(c4_2, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = T, scale ="row",
         show_colnames = T,
         color = colorRampPalette(c("#FFFFFF","#5050f9"))(100))
dev.off()

c4_3 <- read.csv("./4cell_3/learnmodel_inte/mm10_12_chrhmm_overlap.txt", sep="\t")
c4_3 <- c4_3[,2:8]
c4_3 <- c4_3[c(11,1,12,10,9,3,2,4,6,5,7,8),]##与figure2c中的states顺序对应好，图可直接用
c4_3 <- c4_3[,-1]
library(pheatmap)
#c2_1 <- apply(c2_1, 2, function(x) log10(x+0.1))
pdf("./c4_3_state_overlap.pdf", width = 4, height = 6)
pheatmap(c4_3, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = T, scale ="row",
         show_colnames = T,
         color = colorRampPalette(c("#FFFFFF","#5050f9"))(100))
dev.off()

c4_4 <- read.csv("./4cell_4/learnmodel/mm10_12_chrhmm_overlap.txt", sep="\t")
c4_4 <- c4_4[,2:8]
c4_4 <- c4_4[c(1,3,5,4,2,10,12,6,7,8,9,11),]##与figure2c中的states顺序对应好，图可直接用
c4_4 <- c4_4[,-1]
library(pheatmap)
#c2_1 <- apply(c2_1, 2, function(x) log10(x+0.1))
pdf("./c4_4_state_overlap.pdf", width = 4, height = 6)
pheatmap(c4_4, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = T, scale ="row",
         show_colnames = T,
         color = colorRampPalette(c("#FFFFFF","#5050f9"))(100))
dev.off()

