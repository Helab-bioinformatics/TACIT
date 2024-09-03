setwd("/media/helab/helab-7/LM_all/20240427_IVF_embryo_TACIT/IVF_embryo_TACIT")

library(cisTopic)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(tidyr)
library(dplyr)
library(Seurat)
library(Signac)
library(patchwork)
library(stats)

######################################### H3K4me3 ############################################
pathToBams1 <- c('./early2cell_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.2K.windows.bed"
#regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/merge_50kb.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3', min.cells = 0, min.regions = 0)
e2_k4me3_counts <- as.matrix(cisTopicObject_h3k4me3@count.matrix)

pathToBams1 <- c('./late2cell_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.2K.windows.bed"
#regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/merge_50kb.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3', min.cells = 0, min.regions = 0)
l2_k4me3_counts <- as.matrix(cisTopicObject_h3k4me3@count.matrix)

### filter embryos ####
cell_id <- read.csv("./early2cell_k4me3_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
e2_k4me3_counts <- e2_k4me3_counts[,which(colnames(e2_k4me3_counts) %in% cell_id$id)]

cell_id <- read.csv("./late_k4me3_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
l2_k4me3_counts <- l2_k4me3_counts[,which(colnames(l2_k4me3_counts) %in% cell_id$id)]

all_k4me3_counts <- cbind(e2_k4me3_counts, l2_k4me3_counts)
all_k4me3 <- CreateSeuratObject(
  counts = all_k4me3_counts,
  assay = 'peaks',
  project = 'all_k4me3',
  min.cells = 1
)
all_k4me3 <- RunTFIDF(all_k4me3)
all_k4me3 <- FindTopFeatures(all_k4me3, min.cutoff = 'q0')
all_k4me3 <- RunSVD(object = all_k4me3, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("./all2cell_h3k4me3_2kbbin_withdim1.pdf")
DepthCor(all_k4me3)
###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
results <- c()
titration <- c(8, 10, 12, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50)
for (i in 1:length(titration)) {
  all_k4me3 <- RunUMAP(object = all_k4me3, reduction = 'lsi', dims = 1:titration[i], n.neighbors = 10L)#或可调试reduction
  all_k4me3 <- FindNeighbors(object = all_k4me3, reduction = 'lsi', dims = 1:titration[i])
  all_k4me3 <- FindClusters(object = all_k4me3, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
  p1 <- DimPlot(object = all_k4me3, label = TRUE, pt.size = 2)
  p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
  print(p1)
  
  label <- rownames(as.data.frame(all_k4me3@active.ident))
  label <- as.data.frame(label)
  label[grep("E2",label[,1]),2] <- "early2cell"
  label[grep("L2",label[,1]),2] <- "late2cell"
  all_k4me3@meta.data$group.ident <- as.factor(label$V2)
  p2 <- DimPlot(object = all_k4me3, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
  p2 <- p2 + ggtitle(paste0("titration LSI", titration[i]))
  print(p2)
}
dev.off()

saveRDS(all_k4me3, "./seTACIT_h3k4me3_2kbbin.rds")


################################################# late 2cell ###################################################
all_k4me3 <- readRDS("./seTACIT_h3k4me3_5kbbin.rds")
h3k4me3_late <- subset(all_k4me3, cells = colnames(all_k4me3)[grep("late2cell",all_k4me3$group.ident)])
h3k4me3_late <- RunTFIDF(h3k4me3_late)
h3k4me3_late <- FindTopFeatures(h3k4me3_late, min.cutoff = 'q0')
h3k4me3_late <- RunSVD(object = h3k4me3_late, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("./late2cell_h3k4me3_late_5kbbin_titration.pdf")
DepthCor(h3k4me3_late)
###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
results <- c()
titration <- c(3:50)
for (i in 1:length(titration)) {
  h3k4me3_late <- RunUMAP(object = h3k4me3_late, reduction = 'lsi', dims = 2:titration[i], n.neighbors = 13L)#或可调试reduction
  h3k4me3_late <- FindNeighbors(object = h3k4me3_late, reduction = 'lsi', dims = 2:titration[i])
  h3k4me3_late <- FindClusters(object = h3k4me3_late, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
  p1 <- DimPlot(object = h3k4me3_late, label = TRUE, pt.size = 2) + NoLegend()
  p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
  print(p1)

  label <- rownames(as.data.frame(h3k4me3_late@active.ident))
  label <- as.data.frame(label)
  label[grep("E2",label[,1]),2] <- "early2cell"
  label[grep("L2",label[,1]),2] <- "late2cell"
  label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
  label <- label[,c("label", "x5", "x6", "V2")]
  label$cluster <- h3k4me3_late@meta.data[["seurat_clusters"]]
  a <- label[,c("x5", "cluster")]
  library(dplyr)
  same <- a %>%
    group_by(across(everything())) %>%  # 根据所有列对数据进行分组
    filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
    ungroup() 
  same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
  results[i] <- same_ration 
}
names(results) <- paste0("LSI_", titration)
results
dev.off()
write.table(results, "late2cell_h3k4me3_5kbbin_titration.txt", sep = "\t", quote = F)


pdf("late2cell_h3k4me3_2kbbin_same_diff_embryo_ratio.pdf")
h3k4me3_late <- RunUMAP(object = h3k4me3_late, reduction = 'lsi', dims = 2:44, n.neighbors = 13L)#或可调试reduction
h3k4me3_late <- FindNeighbors(object = h3k4me3_late, reduction = 'lsi', dims = 2:44)
h3k4me3_late <- FindClusters(object = h3k4me3_late, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k4me3_late, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k4me3_late@active.ident))
label <- as.data.frame(label)
label[grep("IVF_6",label[,1]),2] <- "IVF_6"
label[grep("IVF_7",label[,1]),2] <- "IVF_7"
h3k4me3_late@meta.data$tab.ident <- as.factor(label$V2)
DimPlot(object = h3k4me3_late, label = TRUE, group.by = "tab.ident",pt.size = 2)

label <- rownames(as.data.frame(h3k4me3_late@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
label$cluster <- h3k4me3_late@meta.data[["seurat_clusters"]]
umap <- as.data.frame(h3k4me3_late@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_2 > 0),3] <- "c1"
umap[which(umap$UMAP_2 <= 0),3] <- "c2"
h3k4me3_late@meta.data$cluster.ident <- as.factor(umap$V3)
DimPlot(object = h3k4me3_late, label = TRUE, group.by = "cluster.ident",pt.size = 2)

a <- data.frame(embryo=label$x5, cluster=umap$V3)
library(dplyr)
same <- a %>%
  group_by(across(everything())) %>%  # 根据所有列对数据进行分组
  filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
  ungroup() 
same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
same_ration

dat <- data.frame(group=c("diff", "same"), var=c(0.40, 0.60))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

ggplot(dat, aes(x="", y=var, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
#dev.off()


# evaluate signals around zga-related genes
pathToBams1 <- c('./late2cell_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/majorZGA.genebody.5kb.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3', min.cells = 0, min.regions = 0)
h3k4me3_zga_late <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
h3k4me3_zga_late <- colSums(h3k4me3_zga_late)

cell_id <- read.csv("./late_k4me3_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k4me3_zga_late <- h3k4me3_zga_late[which(names(h3k4me3_zga_late) %in% cell_id$id)]
h3k4me3_late@meta.data$zga.signal <- as.numeric(log10(h3k4me3_zga_late))

h3k4me3_late_counts <- as.matrix(h3k4me3_late@assays[["peaks"]]@counts)
h3k4me3_late_counts <- colSums(h3k4me3_late_counts)
h3k4me3_zga_late_nor <- h3k4me3_zga_late/h3k4me3_late_counts
h3k4me3_late@meta.data$zga.signal.nor <- as.numeric(h3k4me3_zga_late_nor)

library(ggpubr)
p1 <- VlnPlot(h3k4me3_late, features = "zga.signal", pt.size = 0.1, group.by = "cluster.ident")
p1 <- p1+ stat_compare_means()
p1
p2 <- VlnPlot(h3k4me3_late, features = "zga.signal.nor", pt.size = 0.1, group.by = "cluster.ident")
p2 <- p2+ stat_compare_means()
p2
dev.off()

# find differential peaks between clusters
da_peaks <- FindMarkers(
  object = h3k4me3_late,
  ident.1 = "0",
  ident.2 = "1",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
da_peaks$loci <- rownames(da_peaks)
da_peaks <- separate(da_peaks, loci, c("min", "hh"), sep = ":")
da_peaks <- separate(da_peaks, hh, c("hh", "hhh"), sep = "-")
da_peaks_cluster0 <- da_peaks[which(da_peaks$pct.1 > da_peaks$pct.2 & da_peaks$p_val < 0.01),]
da_peaks_cluster1 <- da_peaks[which(da_peaks$pct.1 < da_peaks$pct.2 & da_peaks$p_val < 0.05),]
write.table(da_peaks_cluster0[,c("min", "hh", "hhh")], "late2cell_h3k4me3_2kbbin_da_peaks_cluster0.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(da_peaks_cluster1[,c("min", "hh", "hhh")], "late2cell_h3k4me3_2kbbin_da_peaks_cluster1.bed", sep = "\t", quote = F, col.names = F, row.names = F)


##################### calculate cell-to-cell distance ######
h3k4me3_late <- readRDS("./h3k4me3_late_2kbbin.rds")
label <- rownames(as.data.frame(h3k4me3_late@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
umap <- as.data.frame(h3k4me3_late@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_2 > 0),3] <- "c1"
umap[which(umap$UMAP_2 <= 0),3] <- "c2"
label$cluster <- as.factor(umap$V3)

h3k4me3_late@meta.data$cluster.ident <- as.factor(umap$V3)
lsi <- as.data.frame(h3k4me3_late@reductions[["lsi"]]@cell.embeddings)
label <- cbind(label, lsi)

embryos <- unique(label$x5) ## 根据embryo的编号提取出同一个embryo的两个细胞，然后计算lsi上的相似性
intra_dis <- c()
for (i in 1:length(embryos)) {
  dat_select <- label[grep(embryos[i], label$x5),]
  lsi_1 <- as.numeric(dat_select[1,6:21])
  lsi_2 <- as.numeric(dat_select[2,6:21])
  cor <- cor(lsi_1, lsi_2, method = "spearman")
  dis <- 1-cor
  intra_dis <- c(dis, intra_dis)
}
intra_dis
intra_dis <- data.frame(order=1:length(intra_dis), intra_dis=sort(intra_dis, decreasing = T))
pdf("h3k4me3_late2cell_intra_embryo_distance.pdf", width = 5, height = 6)
DimPlot(object = h3k4me3_late, label = TRUE, group.by = "cluster.ident",pt.size = 2)
ggplot(intra_dis,aes(x=order,y=intra_dis))+
  geom_point()+
  theme_bw()
dev.off()

dif_embryo <- intra_dis[intra_dis$intra_dis > 0.5, ]
dif_ratio <- nrow(dif_embryo)/length(embryos) # 60%
write.table(label, "h3k4me3_late_label.txt", sep = "\t", quote = F)
write.table(intra_dis, "h3k4me3_late_intra_embryo_distance.txt", sep = "\t", quote = F)

intra_dis_k4_late <- intra_dis

########################################################### early2cell ########################################
h3k4me3_early <- subset(all_k4me3, cells = colnames(all_k4me3)[grep("early2cell",all_k4me3$group.ident)])
h3k4me3_early <- RunTFIDF(h3k4me3_early)
h3k4me3_early <- FindTopFeatures(h3k4me3_early, min.cutoff = 'q0')
h3k4me3_early <- RunSVD(object = h3k4me3_early, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("./early2cell_h3k4me3_early_2kbbin_titration.pdf")
DepthCor(h3k4me3_early)
###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
results <- c()
titration <- c(8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50)
for (i in 1:length(titration)) {
  h3k4me3_early <- RunUMAP(object = h3k4me3_early, reduction = 'lsi', dims = 2:titration[i], n.neighbors = 10L)#或可调试reduction
  h3k4me3_early <- FindNeighbors(object = h3k4me3_early, reduction = 'lsi', dims = 2:titration[i])
  h3k4me3_early <- FindClusters(object = h3k4me3_early, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
  p1 <- DimPlot(object = h3k4me3_early, label = TRUE, pt.size = 2) + NoLegend()
  p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
  print(p1)
  
  label <- rownames(as.data.frame(h3k4me3_early@active.ident))
  label <- as.data.frame(label)
  label[grep("E2",label[,1]),2] <- "early2cell"
  label[grep("L2",label[,1]),2] <- "early2cell"
  label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
  label <- label[,c("label", "x5", "x6", "V2")]
  label$cluster <- h3k4me3_early@meta.data[["seurat_clusters"]]
  a <- label[,c("x5", "cluster")]
  library(dplyr)
  same <- a %>%
    group_by(across(everything())) %>%  # 根据所有列对数据进行分组
    filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
    ungroup() 
  same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
  results[i] <- same_ration 
}
names(results) <- paste0("LSI_", titration)
results
dev.off()
write.table(results, "early2cell_h3k4me3_2kbbin_titration.txt", sep = "\t", quote = F)


pdf("early2cell_h3k4me3_2kbbin_same_diff_embryo_ratio.pdf")
h3k4me3_early <- RunUMAP(object = h3k4me3_early, reduction = 'lsi', dims = 2:18, n.neighbors = 10L)#或可调试reduction
h3k4me3_early <- FindNeighbors(object = h3k4me3_early, reduction = 'lsi', dims = 2:18)
h3k4me3_early <- FindClusters(object = h3k4me3_early, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k4me3_early, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k4me3_early@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
label$cluster <- h3k4me3_early@meta.data[["seurat_clusters"]]
umap <- as.data.frame(h3k4me3_early@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_1 > 0),3] <- "c1"
umap[which(umap$UMAP_1 <= 0),3] <- "c2"

a <- data.frame(embryo=label$x5, cluster=umap$V3)
library(dplyr)
same <- a %>%
  group_by(across(everything())) %>%  # 根据所有列对数据进行分组
  filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
  ungroup() 
same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
same_ration

dat <- data.frame(group=c("diff", "same"), var=c(0.45, 0.55))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

h3k4me3_early@meta.data$cluster.ident <- as.factor(umap$V3)
DimPlot(object = h3k4me3_early, label = TRUE, pt.size = 2, group.by = "cluster.ident") 
ggplot(dat, aes(x="", y=var, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
#dev.off()


# evaluate signals around zga-related genes
pathToBams1 <- c('./early2cell_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/majorZGA.genebody.5kb.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3', min.cells = 0, min.regions = 0)
h3k4me3_zga_early <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
h3k4me3_zga_early <- colSums(h3k4me3_zga_early)

cell_id <- read.csv("./early2cell_k4me3_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k4me3_zga_early <- h3k4me3_zga_early[which(names(h3k4me3_zga_early) %in% cell_id$id)]
h3k4me3_early@meta.data$zga.signal <- as.numeric(log10(h3k4me3_zga_early))

h3k4me3_early_counts <- as.matrix(h3k4me3_early@assays[["peaks"]]@counts)
h3k4me3_early_counts <- colSums(h3k4me3_early_counts)
h3k4me3_zga_early_nor <- h3k4me3_zga_early/h3k4me3_early_counts
h3k4me3_early@meta.data$zga.signal.nor <- as.numeric(h3k4me3_zga_early_nor)

library(ggpubr)
p1 <- VlnPlot(h3k4me3_early, features = "zga.signal", pt.size = 0.1, group.by = "cluster.ident")
p1 <- p1+ stat_compare_means()
p1
p2 <- VlnPlot(h3k4me3_early, features = "zga.signal.nor", pt.size = 0.1, group.by = "cluster.ident")
p2 <- p2+ stat_compare_means()
p2
dev.off()

# find differential peaks between clusters
da_peaks <- FindMarkers(
  object = h3k4me3_early,
  ident.1 = "0",
  ident.2 = "1",
  test.use = 'LR',
  earlynt.vars = 'nCount_peaks'
)
da_peaks$loci <- rownames(da_peaks)
da_peaks <- separate(da_peaks, loci, c("min", "hh"), sep = ":")
da_peaks <- separate(da_peaks, hh, c("hh", "hhh"), sep = "-")
da_peaks_cluster0 <- da_peaks[which(da_peaks$pct.1 > da_peaks$pct.2 & da_peaks$p_val < 0.05),]
da_peaks_cluster1 <- da_peaks[which(da_peaks$pct.1 < da_peaks$pct.2 & da_peaks$p_val < 0.05),]
write.table(da_peaks_cluster0[,c("min", "hh", "hhh")], "early2cell_h3k4me3_2kbbin_da_peaks_cluster0.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(da_peaks_cluster1[,c("min", "hh", "hhh")], "early2cell_h3k4me3_2kbbin_da_peaks_cluster1.bed", sep = "\t", quote = F, col.names = F, row.names = F)


##################### calculate cell-to-cell distance ######
h3k4me3_early <- readRDS("./h3k4me3_early_2kbbin.rds")
label <- rownames(as.data.frame(h3k4me3_early@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
umap <- as.data.frame(h3k4me3_early@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_1 > 0),3] <- "c1"
umap[which(umap$UMAP_1 <= 0),3] <- "c2"
label$cluster <- as.factor(umap$V3)
h3k4me3_early@meta.data$cluster.ident <- as.factor(umap$V3)
lsi <- as.data.frame(h3k4me3_early@reductions[["lsi"]]@cell.embeddings)
label <- cbind(label, lsi)

embryos <- unique(label$x5) ## 根据embryo的编号提取出同一个embryo的两个细胞，然后计算lsi上的相似性
intra_dis <- c()
for (i in 1:length(embryos)) {
  dat_select <- label[grep(embryos[i], label$x5),]
  lsi_1 <- as.numeric(dat_select[1,6:21])
  lsi_2 <- as.numeric(dat_select[2,6:21])
  cor <- cor(lsi_1, lsi_2, method = "spearman")
  dis <- 1-cor
  intra_dis <- c(dis, intra_dis)
}
intra_dis
intra_dis <- data.frame(order=1:length(intra_dis), intra_dis=sort(intra_dis, decreasing = T))
pdf("h3k4me3_early2cell_intra_embryo_distance.pdf", width = 5, height = 6)
DimPlot(object = h3k4me3_early, label = TRUE, pt.size = 2, group.by = "cluster.ident") 
ggplot(intra_dis,aes(x=order,y=intra_dis))+
  geom_point()+
  theme_bw()
dev.off()
write.table(label, "h3k4me3_early_label.txt", sep = "\t", quote = F)
write.table(intra_dis, "h3k4me3_early_intra_embryo_distance.txt", sep = "\t", quote = F)

dif_embryo <- intra_dis[intra_dis$intra_dis > 0.5, ]
dif_ratio <- nrow(dif_embryo)/length(embryos) # 57%
intra_dis_k4_early <- intra_dis

saveRDS(h3k4me3_early, "h3k4me3_early_2kbbin.rds")
saveRDS(h3k4me3_late, "h3k4me3_late_2kbbin.rds")

##############################################################################################
######################################### h3k27ac ############################################
pathToBams1 <- c('./early_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.10K.windows.bed"
#regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/merge_50kb.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
e2_k27ac_counts <- as.matrix(cisTopicObject_h3k27ac@count.matrix)

pathToBams1 <- c('./late2cell_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.10K.windows.bed"
#regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/merge_50kb.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
l2_k27ac_counts <- as.matrix(cisTopicObject_h3k27ac@count.matrix)

### filter embryos ####
cell_id <- read.csv("./early_k27ac_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
e2_k27ac_counts <- e2_k27ac_counts[,which(colnames(e2_k27ac_counts) %in% cell_id$id)]

cell_id <- read.csv("./late_k27acall_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
l2_k27ac_counts <- l2_k27ac_counts[,which(colnames(l2_k27ac_counts) %in% cell_id$id)]

all_k27ac_counts <- cbind(e2_k27ac_counts, l2_k27ac_counts)
all_k27ac <- CreateSeuratObject(
  counts = all_k27ac_counts,
  assay = 'peaks',
  project = 'all_k27ac',
  min.cells = 5
)
all_k27ac <- RunTFIDF(all_k27ac)
all_k27ac <- FindTopFeatures(all_k27ac, min.cutoff = 'q0')
all_k27ac <- RunSVD(object = all_k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("./all2cell_h3k27ac_10kbbin.pdf")
DepthCor(all_k27ac)
###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
results <- c()
titration <- c(8, 10, 12, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50)
for (i in 1:length(titration)) {
  all_k27ac <- RunUMAP(object = all_k27ac, reduction = 'lsi', dims = 2:titration[i], n.neighbors = 10L)#或可调试reduction
  all_k27ac <- FindNeighbors(object = all_k27ac, reduction = 'lsi', dims = 2:titration[i])
  all_k27ac <- FindClusters(object = all_k27ac, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
  p1 <- DimPlot(object = all_k27ac, label = TRUE, pt.size = 2)
  p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
  print(p1)
  
  label <- rownames(as.data.frame(all_k27ac@active.ident))
  label <- as.data.frame(label)
  label[grep("E2",label[,1]),2] <- "early2cell"
  label[grep("L2",label[,1]),2] <- "late2cell"
  all_k27ac@meta.data$group.ident <- as.factor(label$V2)
  p2 <- DimPlot(object = all_k27ac, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
  p2 <- p2 + ggtitle(paste0("titration LSI", titration[i]))
  print(p2)
}
dev.off()

saveRDS(all_k27ac, "./seTACIT_h3k27ac_10kbbin.rds")



################################################# late 2cell ###################################################
set.seed(12)
all_k27ac <- readRDS("./seTACIT_h3k27ac_5kbin.rds")
h3k27ac_late <- subset(all_k27ac, cells = colnames(all_k27ac)[grep("late2cell",all_k27ac$group.ident)])
h3k27ac_late <- RunTFIDF(h3k27ac_late)
h3k27ac_late <- FindTopFeatures(h3k27ac_late, min.cutoff = 'q0')
h3k27ac_late <- RunSVD(object = h3k27ac_late, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("./late2cell_h3k27ac_late_5kbbin_titration.pdf")
DepthCor(h3k27ac_late)
###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
results <- c()
titration <- c(5:50)
for (i in 1:length(titration)) {
  h3k27ac_late <- RunUMAP(object = h3k27ac_late, reduction = 'lsi', dims = 2:titration[i], n.neighbors = 15L)#或可调试reduction
  h3k27ac_late <- FindNeighbors(object = h3k27ac_late, reduction = 'lsi', dims = 2:titration[i])
  h3k27ac_late <- FindClusters(object = h3k27ac_late, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
  p1 <- DimPlot(object = h3k27ac_late, label = TRUE, pt.size = 2) + NoLegend()
  p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
  print(p1)
  
  label <- rownames(as.data.frame(h3k27ac_late@active.ident))
  label <- as.data.frame(label)
  label[grep("E2",label[,1]),2] <- "early2cell"
  label[grep("L2",label[,1]),2] <- "late2cell"
  label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
  label <- label[,c("label", "x5", "x6", "V2")]
  label$cluster <- h3k27ac_late@meta.data[["seurat_clusters"]]
  a <- label[,c("x5", "cluster")]
  library(dplyr)
  same <- a %>%
    group_by(across(everything())) %>%  # 根据所有列对数据进行分组
    filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
    ungroup() 
  same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
  results[i] <- same_ration 
}
names(results) <- paste0("LSI_", titration)
dev.off()
results
write.table(results, "late2cell_h3k27ac_late_5kbbin_titration.txt", sep = "\t", quote = F)


pdf("late2cell_h3k27ac_5kbbin_late_same_diff_embryo_ratio.pdf")
h3k27ac_late <- RunUMAP(object = h3k27ac_late, reduction = 'lsi', dims = 2:12, n.neighbors = 10L)#或可调试reduction
h3k27ac_late <- FindNeighbors(object = h3k27ac_late, reduction = 'lsi', dims = 2:12)
h3k27ac_late <- FindClusters(object = h3k27ac_late, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k27ac_late, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k27ac_late@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "late2cell"
label[grep("L2",label[,1]),2] <- "late2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
label$cluster <- h3k27ac_late@meta.data[["seurat_clusters"]]
umap <- as.data.frame(h3k27ac_late@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_1 > -2),3] <- "c1"
umap[which(umap$UMAP_1 <= -2),3] <- "c2"
a <- data.frame(embryo=label$x5, cluster=label$cluster)
library(dplyr)
same <- a %>%
  group_by(across(everything())) %>%  # 根据所有列对数据进行分组
  filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
  ungroup() 
same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
same_ration

dat <- data.frame(group=c("diff", "same"), var=c(0.41, 0.59))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

DimPlot(object = h3k27ac_late, label = TRUE, pt.size = 2) + NoLegend()
ggplot(dat, aes(x="", y=var, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
#dev.off()

# evaluate signals around zga-related genes
pathToBams1 <- c('./late2cell_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/majorZGA.genebody.5kb.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='late2cell_h3k27ac', min.cells = 0, min.regions = 0)
h3k27ac_zga_late <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
h3k27ac_zga_late <- colSums(h3k27ac_zga_late)

cell_id <- read.csv("./late_k27acall_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k27ac_zga_late <- h3k27ac_zga_late[which(names(h3k27ac_zga_late) %in% cell_id$id)]
h3k27ac_late@meta.data$zga.signal <- as.numeric(log10(h3k27ac_zga_late))

h3k27ac_late_counts <- as.matrix(h3k27ac_late@assays[["peaks"]]@counts)
h3k27ac_late_counts <- colSums(h3k27ac_late_counts)
h3k27ac_zga_late_nor <- h3k27ac_zga_late/h3k27ac_late_counts
h3k27ac_late@meta.data$zga.signal.nor <- as.numeric(h3k27ac_zga_late_nor)

library(ggpubr)
p1 <- VlnPlot(h3k27ac_late, features = "zga.signal", pt.size = 0.1)
p1 <- p1+ stat_compare_means()
p1
p2 <- VlnPlot(h3k27ac_late, features = "zga.signal.nor", pt.size = 0.1)
p2 <- p2+ stat_compare_means()
p2
dev.off()

# find differential peaks between clusters
da_peaks <- FindMarkers(
  object = h3k27ac_late,
  ident.1 = "0",
  ident.2 = "1",
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
da_peaks$loci <- rownames(da_peaks)
da_peaks <- separate(da_peaks, loci, c("min", "hh"), sep = ":")
da_peaks <- separate(da_peaks, hh, c("hh", "hhh"), sep = "-")
da_peaks_cluster0 <- da_peaks[which(da_peaks$pct.1 > da_peaks$pct.2 & da_peaks$p_val < 0.05),]
da_peaks_cluster1 <- da_peaks[which(da_peaks$pct.1 < da_peaks$pct.2 & da_peaks$p_val < 0.05),]
write.table(da_peaks_cluster0[,c("min", "hh", "hhh")], "late2cell_h3k27ac_5kbbin_da_peaks_cluster0.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(da_peaks_cluster1[,c("min", "hh", "hhh")], "late2cell_h3k27ac_5kbbin_da_peaks_cluster1.bed", sep = "\t", quote = F, col.names = F, row.names = F)


##################### calculate cell-to-cell distance ######
h3k27ac_late <- readRDS("./h3k27ac_late_5kbbin.rds")
label <- rownames(as.data.frame(h3k27ac_late@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
umap <- as.data.frame(h3k27ac_late@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_2 > -2),3] <- "c1"
umap[which(umap$UMAP_2 <= -2),3] <- "c2"
label$cluster <- as.factor(umap$V3)
h3k27ac_late@meta.data$cluster.ident <- as.factor(umap$V3)
lsi <- as.data.frame(h3k27ac_late@reductions[["lsi"]]@cell.embeddings)
label <- cbind(label, lsi)

embryos <- unique(label$x5) ## 根据embryo的编号提取出同一个embryo的两个细胞，然后计算lsi上的相似性
intra_dis <- c()
for (i in 1:length(embryos)) {
  dat_select <- label[grep(embryos[i], label$x5),]
  lsi_1 <- as.numeric(dat_select[1,6:21])
  lsi_2 <- as.numeric(dat_select[2,6:21])
  cor <- cor(lsi_1, lsi_2, method = "spearman")
  dis <- 1-cor
  intra_dis <- c(dis, intra_dis)
}
intra_dis
intra_dis <- data.frame(order=1:length(intra_dis), intra_dis=sort(intra_dis, decreasing = T))
pdf("h3k27ac_late2cell_intra_embryo_distance.pdf", width = 5, height = 6)
DimPlot(object = h3k27ac_late, label = TRUE, pt.size = 2, group.by = "cluster.ident")
ggplot(intra_dis,aes(x=order,y=intra_dis))+
  geom_point()+
  theme_bw()
dev.off()

dif_embryo <- intra_dis[intra_dis$intra_dis > 0.5, ]
dif_ratio <- nrow(dif_embryo)/length(embryos) # 67%
write.table(label, "h3k27ac_late_label.txt", sep = "\t", quote = F)
intra_dis_k27_late <- intra_dis

########################################################### early2cell ########################################
all_k27ac <- readRDS("./seTACIT_h3k27ac_5kbin.rds")
h3k27ac_early <- subset(all_k27ac, cells = colnames(all_k27ac)[grep("early2cell",all_k27ac$group.ident)])
h3k27ac_early <- RunTFIDF(h3k27ac_early)
h3k27ac_early <- FindTopFeatures(h3k27ac_early, min.cutoff = 'q0')
h3k27ac_early <- RunSVD(object = h3k27ac_early, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')

pdf("./early2cell_h3k27ac_early_5kbbin_titration.pdf")
DepthCor(h3k27ac_early)
###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
results <- c()
titration <- c(3:50)
for (i in 1:length(titration)) {
  h3k27ac_early <- RunUMAP(object = h3k27ac_early, reduction = 'lsi', dims = 2:titration[i], n.neighbors = 13L)#或可调试reduction
  h3k27ac_early <- FindNeighbors(object = h3k27ac_early, reduction = 'lsi', dims = 2:titration[i])
  h3k27ac_early <- FindClusters(object = h3k27ac_early, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
  p1 <- DimPlot(object = h3k27ac_early, label = TRUE, pt.size = 2) + NoLegend()
  p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
  print(p1)
  
  label <- rownames(as.data.frame(h3k27ac_early@active.ident))
  label <- as.data.frame(label)
  label[grep("E2",label[,1]),2] <- "early2cell"
  label[grep("L2",label[,1]),2] <- "early2cell"
  label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
  label <- label[,c("label", "x5", "x6", "V2")]
  label$cluster <- h3k27ac_early@meta.data[["seurat_clusters"]]
  a <- label[,c("x5", "cluster")]
  library(dplyr)
  same <- a %>%
    group_by(across(everything())) %>%  # 根据所有列对数据进行分组
    filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
    ungroup() 
  same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
  results[i] <- same_ration 
}
names(results) <- paste0("LSI_", titration)
results
dev.off()
write.table(results, "early2cell_h3k27ac_5kbbin_titration.txt", sep = "\t", quote = F)


pdf("early2cell_h3k27ac_5kbbin_same_diff_embryo_ratio.pdf")
h3k27ac_early <- RunUMAP(object = h3k27ac_early, reduction = 'lsi', dims = 2:13, n.neighbors = 13L)#或可调试reduction
h3k27ac_early <- FindNeighbors(object = h3k27ac_early, reduction = 'lsi', dims = 2:13)
h3k27ac_early <- FindClusters(object = h3k27ac_early, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k27ac_early, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k27ac_early@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
label$cluster <- h3k27ac_early@meta.data[["seurat_clusters"]]
umap <- as.data.frame(h3k27ac_early@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_1 > 2),3] <- "c1"
umap[which(umap$UMAP_1 <= 2),3] <- "c2"
a <- data.frame(embryo=label$x5, cluster=umap$V3)
library(dplyr)
same <- a %>%
  group_by(across(everything())) %>%  # 根据所有列对数据进行分组
  filter(n() > 1) %>%  # 仅保留在组内出现多于一次的行
  ungroup() 
same_ration <- nrow(same)/nrow(a)  # 数据框的总行数
same_ration

label <- rownames(as.data.frame(h3k27ac_early@active.ident))
label <- as.data.frame(label)
label[grep("IVF_1",label[,1]),2] <- "IVF_1"
label[grep("IVF_2",label[,1]),2] <- "IVF_2"
h3k27ac_early@meta.data$tab.ident <- as.factor(label$V2)
DimPlot(object = h3k27ac_early, label = TRUE, group.by = "tab.ident",pt.size = 2)

h3k27ac_early@meta.data$cluster.ident <- as.factor(umap$V3)
DimPlot(object = h3k27ac_early, label = TRUE, group.by = "cluster.ident",pt.size = 2)

dat <- data.frame(group=c("diff", "same"), var=c(0.43, 0.57))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

ggplot(dat, aes(x="", y=var, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
#dev.off()


# evaluate signals around zga-related genes
pathToBams1 <- c('./early_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/majorZGA.genebody.5kb.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
h3k27ac_zga_early <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
h3k27ac_zga_early <- colSums(h3k27ac_zga_early)

cell_id <- read.csv("./early_k27ac_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k27ac_zga_early <- h3k27ac_zga_early[which(names(h3k27ac_zga_early) %in% cell_id$id)]
h3k27ac_early@meta.data$zga.signal <- as.numeric(log10(h3k27ac_zga_early))

h3k27ac_early_counts <- as.matrix(h3k27ac_early@assays[["peaks"]]@counts)
h3k27ac_early_counts <- colSums(h3k27ac_early_counts)
h3k27ac_zga_early_nor <- h3k27ac_zga_early/h3k27ac_early_counts
h3k27ac_early@meta.data$zga.signal.nor <- as.numeric(h3k27ac_zga_early_nor)

library(ggpubr)
p1 <- VlnPlot(h3k27ac_early, features = "zga.signal", pt.size = 0.1, group.by = "cluster.ident")
p1 <- p1+ stat_compare_means()
p1
p2 <- VlnPlot(h3k27ac_early, features = "zga.signal.nor", pt.size = 0.1, group.by = "cluster.ident",)
p2 <- p2+ stat_compare_means()
p2
dev.off()

# find differential peaks between clusters
da_peaks <- FindMarkers(
  object = h3k27ac_early,
  ident.1 = "0",
  ident.2 = "1",
  test.use = 'LR',
  earlynt.vars = 'nCount_peaks'
)
da_peaks$loci <- rownames(da_peaks)
da_peaks <- separate(da_peaks, loci, c("min", "hh"), sep = ":")
da_peaks <- separate(da_peaks, hh, c("hh", "hhh"), sep = "-")
da_peaks_cluster0 <- da_peaks[which(da_peaks$pct.1 > da_peaks$pct.2 & da_peaks$p_val < 0.05),]
da_peaks_cluster1 <- da_peaks[which(da_peaks$pct.1 < da_peaks$pct.2 & da_peaks$p_val < 0.05),]
write.table(da_peaks_cluster0[,c("min", "hh", "hhh")], "early2cell_h3k27ac_2kbbin_da_peaks_cluster0.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(da_peaks_cluster1[,c("min", "hh", "hhh")], "early2cell_h3k27ac_2kbbin_da_peaks_cluster1.bed", sep = "\t", quote = F, col.names = F, row.names = F)

saveRDS(h3k27ac_early, "h3k27ac_early_5kbbin.rds")
saveRDS(h3k27ac_late, "h3k27ac_late_5kbbin.rds")

##################### calcuearly cell-to-cell distance ######
h3k27ac_early <- readRDS("./h3k27ac_early_5kbbin.rds")
label <- rownames(as.data.frame(h3k27ac_early@active.ident))
label <- as.data.frame(label)
label[grep("E2",label[,1]),2] <- "early2cell"
label[grep("L2",label[,1]),2] <- "early2cell"
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6", "V2")]
umap <- as.data.frame(h3k27ac_early@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_2 > 2),3] <- "c1"
umap[which(umap$UMAP_2 <= 2),3] <- "c2"
label$cluster <- as.factor(umap$V3)
h3k27ac_early@meta.data$cluster.ident <- as.factor(umap$V3)
lsi <- as.data.frame(h3k27ac_early@reductions[["lsi"]]@cell.embeddings)
label <- cbind(label, lsi)

embryos <- unique(label$x5) ## 根据embryo的编号提取出同一个embryo的两个细胞，然后计算lsi上的相似性
intra_dis <- c()
for (i in 1:length(embryos)) {
  dat_select <- label[grep(embryos[i], label$x5),]
  lsi_1 <- as.numeric(dat_select[1,6:21])
  lsi_2 <- as.numeric(dat_select[2,6:21])
  cor <- cor(lsi_1, lsi_2, method = "spearman")
  dis <- 1-cor
  intra_dis <- c(dis, intra_dis)
}
intra_dis
intra_dis <- data.frame(order=1:length(intra_dis), intra_dis=sort(intra_dis, decreasing = T))
pdf("h3k27ac_early2cell_intra_embryo_distance.pdf", width = 5, height = 6)
DimPlot(object = h3k27ac_early, label = TRUE, pt.size = 2, group.by = "cluster.ident")
ggplot(intra_dis,aes(x=order,y=intra_dis))+
  geom_point()+
  theme_bw()
dev.off()

dif_embryo <- intra_dis[intra_dis$intra_dis > 0.5, ]
dif_ratio <- nrow(dif_embryo)/length(embryos) # 52%
write.table(label, "h3k27ac_early_label.txt", sep = "\t", quote = F)

intra_dis_k27_early <- intra_dis

##### plot all distance ####
all_distance <- data.frame(distance=c(intra_dis_k4_early$intra_dis, intra_dis_k4_late$intra_dis,
                                      intra_dis_k27_early$intra_dis, intra_dis_k27_late$intra_dis),
                           group=c(rep("h3k4me3_early", nrow(intra_dis_k4_early)), 
                                   rep("h3k4me3_late", nrow(intra_dis_k4_late)),
                                   rep("h3k27ac_early", nrow(intra_dis_k27_early)),
                                   rep("h3k27ac_late", nrow(intra_dis_k27_late))))

pdf("IVF_all_intra_embryo_distance.pdf", width = 8)
ggplot(data=all_distance,mapping = aes(x = group, y = distance))+
  geom_violin(aes(fill = group), trim = T) + 
  geom_boxplot(width = 0.1)+
 # geom_point(position = "jitter")
  theme(legend.position = "none")+
  theme_bw()
dev.off()
write.table(all_distance, "IVF_all_intra_embryo_distance.txt", sep = "\t", quote = F)


####################### plot with in vivo data #######################
dis_ivf <- read.csv("./IVF_all_intra_embryo_distance.txt", sep = "\t")
dis_in_vivo <- read.csv("/media/helab/data3/min/TACIT/20240103_mouse_embryo_TACIT/all_intra_embryo_distance.txt", sep = "\t")
dis_ivf$stage <- "ivf"
dis_in_vivo$stage <- "Invivo"
all <- rbind(dis_in_vivo, dis_ivf)
all$label <- paste(all$stage, all$group, sep = "_")
#all$label <- factor(all$label, levels = c())
pdf("IVF_and_in_vivo_intra_embryo_distance.pdf", width = 8)
ggplot(data=all,mapping = aes(x = label, y = distance))+
  geom_violin(aes(fill = label), trim = T) + 
  geom_boxplot(width = 0.1)+
  # geom_point(position = "jitter")
  theme(legend.position = "none")+
  theme_bw()
dev.off()


############################### distance for mESCs ##################################
pathToBams1 <- c('/media/helab/data1/min/02_tacit/01_cell_lines/20220405_mESC_k4k27k36/1.rawdata/ES_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.10K.windows.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='late2cell_h3k4me3', min.cells = 0, min.regions = 0)
h3k4me3_ES <- as.matrix(cisTopicObject_h3k4me3@count.matrix)

h3k4me3_ES <- CreateSeuratObject(
  counts = h3k4me3_ES,
  assay = 'peaks',
  project = 'all_k27ac',
  min.cells = 5
)
h3k4me3_ES <- RunTFIDF(h3k4me3_ES)
h3k4me3_ES <- FindTopFeatures(h3k4me3_ES, min.cutoff = 'q0')
h3k4me3_ES <- RunSVD(object = h3k4me3_ES, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
h3k4me3_ES_lsi <- as.data.frame(h3k4me3_ES@reductions[["lsi"]]@cell.embeddings)
h3k4me3_ES_lsi <- h3k4me3_ES_lsi[,2:20]

##### spearman cor between all h3k4me3_ES data #########
h3k4me3_ES_spear <- vector("list", nrow(h3k4me3_ES_lsi))
for (i in 1:nrow(h3k4me3_ES_lsi)) {
  row_h3k4me3_ES_lsi <- as.numeric(h3k4me3_ES_lsi[i, ])
  spearmans <- numeric(nrow(h3k4me3_ES_lsi))
  for (j in 1:nrow(h3k4me3_ES_lsi)) {
    row_h3k4me3_ES_lsi_2 <- as.numeric(h3k4me3_ES_lsi[j, ])
    spearmans[j] <- cor(row_h3k4me3_ES_lsi, row_h3k4me3_ES_lsi_2, method = "spearman")
  }
  h3k4me3_ES_spear[[i]] <- spearmans
}
h3k4me3_ES_spear_df <- do.call(rbind, h3k4me3_ES_spear)
h3k4me3_ES_spear_df <- as.numeric(h3k4me3_ES_spear_df)
h3k4me3_ES_spear_df <- h3k4me3_ES_spear_df[h3k4me3_ES_spear_df <1]
h3k4me3_ES_dis <- 1-h3k4me3_ES_spear_df

# ##############################################################################################
# ######################################### H3k27ac ############################################
# pathToBams1 <- c('./early2cell_h3k27ac/')
# bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
# bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
# regions <- "/media/helab/data1/min/00_reference/mm10.2K.windows.bed"
# cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
# e2_k27ac_counts <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
# 
# pathToBams1 <- c('./late2cell_h3k27ac/')
# bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
# bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
# regions <- "/media/helab/data1/min/00_reference/mm10.2K.windows.bed"
# cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
# l2_k27ac_counts <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
# 
# ### filter embryos ####
# cell_id <- read.csv("./early2cell_k27ac_all_QC.txt", sep = "\t", header = F)
# cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
# e2_k27ac_counts <- e2_k27ac_counts[,which(colnames(e2_k27ac_counts) %in% cell_id$id)]
# 
# cell_id <- read.csv("./late_k27ac_all_QC.txt", sep = "\t", header = F)
# cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
# l2_k27ac_counts <- l2_k27ac_counts[,which(colnames(l2_k27ac_counts) %in% cell_id$id)]
# 
# all_k27ac_counts <- cbind(e2_k27ac_counts, l2_k27ac_counts)
# all_k27ac <- CreateSeuratObject(
#   counts = all_k27ac_counts,
#   assay = 'peaks',
#   project = 'all_k27ac',
#   min.cells = 1
# )
# all_k27ac <- RunTFIDF(all_k27ac)
# all_k27ac <- FindTopFeatures(all_k27ac, min.cutoff = 'q0')
# all_k27ac <- RunSVD(object = all_k27ac, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
# 
# pdf("./all2cell_h3k27ac_2kbbin.pdf")
# DepthCor(all_k27ac)
# ###### 通过titration调整胚胎分布在不同cluster的比例，尽可能多一些####
# results <- c()
# titration <- c(8, 10, 12, 14, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50)
# for (i in 1:length(titration)) {
#   all_k27ac <- RunUMAP(object = all_k27ac, reduction = 'lsi', dims = 2:titration[i], n.neighbors = 10L)#或可调试reduction
#   all_k27ac <- FindNeighbors(object = all_k27ac, reduction = 'lsi', dims = 2:titration[i])
#   all_k27ac <- FindClusters(object = all_k27ac, algorithm = 3, resolution = 0.3, verbose = FALSE)#或可调试reduction和algorithm
#   p1 <- DimPlot(object = all_k27ac, label = TRUE, pt.size = 2, title=paste0("titration LSI", titration[i]))
#   p1 <- p1 + ggtitle(paste0("titration LSI", titration[i]))
#   print(p1)
#   
#   label <- rownames(as.data.frame(all_k27ac@active.ident))
#   label <- as.data.frame(label)
#   label[grep("E2",label[,1]),2] <- "early2cell"
#   label[grep("L2",label[,1]),2] <- "late2cell"
#   all_k27ac@meta.data$group.ident <- as.factor(label$V2)
#   p2 <- DimPlot(object = all_k27ac, group.by = "group.ident" ,label = TRUE, pt.size = 2,  title=paste0("titration LSI", titration[i])) 
#   p2 <- p2 + ggtitle(paste0("titration LSI", titration[i]))
#   print(p2)
# }
# dev.off()
# 
# saveRDS(all_k27ac, "./seTACIT_h3k27ac.rds")





