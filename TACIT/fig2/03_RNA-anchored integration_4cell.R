
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

############################################# RNA ########################################################
RNA <- readRDS("./all_RNA.rds")
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

#######################################################################################################################
####################################### 4cell integration #############################################################
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
out.id <- cutree(out.hclust, k=4)

pdf("./RNA_4cell_pseudo.hcluster.pdf")
plot(out.hclust)
dev.off()

############################################# H3K4me3 (as an example) ####################################################
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
h3k4me3 <- RunUMAP(object = h3k4me3, reduction = 'lsi', dims = c(2:15))
h3k4me3 <- FindNeighbors(object = h3k4me3, reduction = 'lsi', dims = c(2:15))
h3k4me3 <- FindClusters(object = h3k4me3, algorithm = 3, resolution = 1, verbose = FALSE)
h3k4me3_count <- as.matrix(h3k4me3@assays[["peaks"]]@counts)
depth <- apply(h3k4me3_count, 2, sum)
h3k4me3$depth <- log10(depth)

pdf("./h3k4me3_4cell.pdf")
DepthCor(h3k4me3)
DimPlot(object = h3k4me3, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(h3k4me3, reduction = "umap", group.by = "group.ident",label = TRUE, pt.size = 2)
FeaturePlot(h3k4me3, reduction = "umap", features = 'depth')
dev.off()

### build gene matrix use signals around TSS regions ####
pathToBams1 <- c('/media/helab/data1/min/02_tacit/03_early_stage/all_h3k4me3/allBam/Bam2000/all/bam/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bai|ArchR|cluster|bed",invert = T, bamFiles1)]
#regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl2kb_all+.bed"
regions <- "/media/helab/data1/min/00_reference/mm10_TSS_fl5kb_all+.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='h3k4me3', min.regions = 0)
h3k4me3_pr <- as.matrix(cisTopicObject@count.matrix)

# use CreatGeneActivityMatrix.R to generate gene activity matrix

geneActivity <- read.csv("./h3k4me3_allstage_tss2kb_ActivityMatrix.txt", sep = " ",check.names=F)
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

####################################### only retain cells with  prediction.score.max > 0.2 across all motalities #######################################################
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
h3k4me1.sh <- paste0("cp ./", rownames(h3k4me1_pre1), " ./4cell/", h3k4me1_pre1$predicted.id, "_inte.bam")
write.table(h3k4me1.sh, "./h3k4me1/h3k4me1_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k4me3_pre1 <- h3k4me3_pre1 [which(h3k4me3_pre1 $predicted.id %in% inte),]
h3k4me3_pre1 <- h3k4me3_pre1[order(h3k4me3_pre1$prediction.score.max, decreasing = T),]
h3k4me3_pre1  <- h3k4me3_pre1 [!duplicated(h3k4me3_pre1 $predicted.id),]
h3k4me3.sh <- paste0("cp ./", rownames(h3k4me3_pre1), " ./4cell/", h3k4me3_pre1$predicted.id, "_inte.bam")
write.table(h3k4me3.sh, "./h3k4me3/h3k4me3_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k27ac_pre1 <- h3k27ac_pre1 [which(h3k27ac_pre1 $predicted.id %in% inte),]
h3k27ac_pre1 <- h3k27ac_pre1[order(h3k27ac_pre1$prediction.score.max, decreasing = T),]
h3k27ac_pre1  <- h3k27ac_pre1 [!duplicated(h3k27ac_pre1 $predicted.id),]
h3k27ac.sh <- paste0("cp ./", rownames(h3k27ac_pre1), " ./4cell/", h3k27ac_pre1$predicted.id, "_inte.bam")
write.table(h3k27ac.sh, "./h3k27ac/h3k27ac_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

h3k36me3_pre1 <- h3k36me3_pre1 [which(h3k36me3_pre1 $predicted.id %in% inte),]
h3k36me3_pre1 <- h3k36me3_pre1[order(h3k36me3_pre1$prediction.score.max, decreasing = T),]
h3k36me3_pre1  <- h3k36me3_pre1 [!duplicated(h3k36me3_pre1 $predicted.id),]
h3k36me3.sh <- paste0("cp ./", rownames(h3k36me3_pre1), " ./4cell/", h3k36me3_pre1$predicted.id, "_inte.bam")
write.table(h3k36me3.sh, "./h3k36me3/h3k36me3_4cell_linked.sh", sep = "\t", quote = F, col.names = F, row.names = F)

cotacit_pre1 <- cotacit_pre1 [which(cotacit_pre1 $predicted.id %in% inte),]
cotacit_pre1 <- cotacit_pre1[order(cotacit_pre1$prediction.score.max, decreasing = T),]
cotacit_pre1  <- cotacit_pre1 [!duplicated(cotacit_pre1 $predicted.id),]
cotacit_pre1$ID <- rownames(cotacit_pre1)
cotacit_pre1 <- separate(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10"), sep = "_")
cotacit_pre1 <- unite(cotacit_pre1, ID, c("v1", "v2", "v3", "v4", "v5", "v6"), sep = "_")
cotacit_pre1$ID2 <- paste0(cotacit_pre1$ID, "_combined_all_h3k27me3.mm10_rmdup.bam")
cotacit_pre1$ID3 <- paste0(cotacit_pre1$ID, "_all_h3k9me3_rmdup.bam")
cotacit.sh <- paste0("cp ./cotacit_data/h3k27me3/", cotacit_pre1$ID2, " ./4cell/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k27me3/cotacit_4cell_linked_h3k27me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
cotacit.sh <- paste0("cp ./cotacit_data/h3k9me3/", cotacit_pre1$ID3, " ./4cell/", cotacit_pre1$predicted.id, "_inte.bam")
write.table(cotacit.sh, "./cotacit_data/integration/h3k9me3/cotacit_4cell_linked_h3k9me3.sh", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(cotacit_pre1[,c("predicted.id","prediction.score.max", "ID2", "ID3")], "./cotacit_data/integration/cotacit_4cell_integration.txt")


############################## annotation chromatin states ########################################
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
