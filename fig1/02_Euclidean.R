rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F) 
################################# calculate Euclidean distance #################################
library(cisTopic)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(Seurat)
library(Signac)


######################### h3k4me3 (as an example) ##########################
setwd("./downsample")
pathToBams1 <- c('./')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|agg|ArchR",invert = T, bamFiles1)]
regions <- "./all_h3k4me3-0.05_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me3')
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(5:20), seed=987, nCores=30, iterations = 500, addModels=FALSE)
pdf("./selectModel.pdf")
cisTopicObject <- selectModel(cisTopicObject)
dev.off()

cell.names <- cisTopicObject@cell.names
label <- as.data.frame(cell.names)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("E2.0|8cell",label[,1]),2] <- "8cell"
label[grep("E2.5|morula",label[,1]),2] <- "morula"
label[grep("E3.0|E3.5|blastocyst",label[,1]),2] <- "blastocyst"
cisTopicObject@cell.data$stage <- label[,2]

#clustering
cisTopicObject <- runtSNE(cisTopicObject, target='cell', seed=123, pca=FALSE, method='Probability')
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, method='Probability')

pdf('./h3k4me3_umap.pdf')
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c('stage'), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, colsVar=(c("#2E94B9","#FFFDC0","#F0B775","#D25565")), intervals=20, labels = T)
dev.off()

h3k4me3_umap <- as.matrix(cisTopicObject@dr[["cell"]][["Umap"]])
umap_2cell <- h3k4me3_umap[grep("2cell", rownames(h3k4me3_umap)),]
umap_2cell <- as.data.frame(umap_2cell)
umap_2cell_out <- umap_2cell[which(umap_2cell$UMAP2 >-5),]

h3k4me3_umap <- as.data.frame(h3k4me3_umap)
h3k4me3_umap <- h3k4me3_umap[!(rownames(h3k4me3_umap) %in% rownames(umap_2cell_out)),]
dis_zygote <- proxy::dist(h3k4me3_umap[grep("zygote", rownames(h3k4me3_umap)),], method = "Euclidean")
dis_2cell <- proxy::dist(h3k4me3_umap[grep("2cell", rownames(h3k4me3_umap)),], method = "Euclidean")
dis_4cell <- proxy::dist(h3k4me3_umap[grep("4cell", rownames(h3k4me3_umap)),], method = "Euclidean")
dis_8cell <- proxy::dist(h3k4me3_umap[grep("8cell", rownames(h3k4me3_umap)),], method = "Euclidean")
dis_morula <- proxy::dist(h3k4me3_umap[grep("morula", rownames(h3k4me3_umap)),], method = "Euclidean")
dis_blastocyst <- proxy::dist(h3k4me3_umap[grep("blastocyst|E3.5|E3.0", rownames(h3k4me3_umap)),], method = "Euclidean")

dis_2cell <- dis_2cell[dis_2cell<15]
max <- median(dis_zygote)
dat <- data.frame(group=c(rep("zygote", length(dis_zygote)), rep("2cell", length(dis_2cell)), rep("4cell", length(dis_4cell)), rep("8cell", length(dis_8cell)), rep("morula", length(dis_morula)), rep("blastocyst", length(dis_blastocyst))),
                  distance=c(dis_zygote/max, dis_2cell/max, dis_4cell/max, dis_8cell/max, dis_morula/max, dis_blastocyst/max))
dat$group <- factor(dat$group, levels = c("zygote", "2cell", "4cell", "8cell", "morula", "blastocyst"))

median <- data.frame(median=c(median(dis_zygote/max), median(dis_2cell/max), median(dis_4cell/max), median(dis_8cell/max), median(dis_morula/max), median(dis_blastocyst/max)),
                     group=c("zygote", "2cell", "4cell", "8cell", "morula", "blastocyst"))

pdf("h3k4me3_euclidean.pdf")
ggplot(dat,mapping = aes(x = group, y = distance))+
  geom_boxplot(aes(fill = group),width = 0.6)+
  theme_bw()+
  geom_text(data=median, aes(x=group, y=median, label=median), nudge_y=1)
dev.off()


################################################# h3k27ac (as an example) #########################################################
setwd("./downsample")
pathToBams1 <- c('./')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("bad|txt|06_reads|bed|cluster|early_2cell",invert = T, bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mm10.10K.windows.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac')

c <- as.matrix(cisTopicObject@count.matrix) 
h3k27ac <- CreateSeuratObject(counts = c,assay = 'peaks',project = 'h3k27ac', min.cells = 5)
h3k27ac <- RunTFIDF(h3k27ac)
h3k27ac <- FindTopFeatures(h3k27ac, min.cutoff = 'q0')
h3k27ac <- RunSVD(object = h3k27ac,assay = 'peaks',reduction.key = 'LSI_',reduction.name = 'lsi')

label <- rownames(as.data.frame(h3k27ac@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("E2.0|8cell",label[,1]),2] <- "8cell"
label[grep("E2.5|morula",label[,1]),2] <- "morula"
label[grep("E3.0|E3.5|blastocyst",label[,1]),2] <- "blastocyst"
h3k27ac@meta.data$group.ident <- as.factor(label[,2])

pdf("./h3k27ac_umap.pdf")
DepthCor(h3k27ac)
h3k27ac <- RunUMAP(object = h3k27ac,reduction = 'lsi',dims = c(2:25))#或可调试reduction
h3k27ac <- FindNeighbors(object = h3k27ac,reduction = 'lsi',dims = c(2:25))
h3k27ac <- FindClusters(object = h3k27ac,algorithm = 3, resolution = 0.3,verbose = FALSE)#或可调试reduction和algorithm
DimPlot(object = h3k27ac, label = TRUE, pt.size = 2) + NoLegend()
DimPlot(object = h3k27ac, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
dev.off()


h3k27ac_umap <- as.data.frame(h3k27ac@reductions$umap@cell.embeddings)
h3k27ac_umap_out <- h3k27ac_umap[which(h3k27ac_umap$UMAP_1 < -4 & h3k27ac_umap$UMAP_2 < -4),]
h3k27ac_umap <- h3k27ac_umap[!(rownames(h3k27ac_umap) %in% rownames(h3k27ac_umap_out)),]
dis_zygote <- proxy::dist(h3k27ac_umap[grep("zygote", rownames(h3k27ac_umap)),], method = "Euclidean")
dis_2cell <- proxy::dist(h3k27ac_umap[grep("2cell", rownames(h3k27ac_umap)),], method = "Euclidean")
dis_4cell <- proxy::dist(h3k27ac_umap[grep("4cell", rownames(h3k27ac_umap)),], method = "Euclidean")
dis_8cell <- proxy::dist(h3k27ac_umap[grep("8cell", rownames(h3k27ac_umap)),], method = "Euclidean")
dis_morula <- proxy::dist(h3k27ac_umap[grep("morula|E2.5", rownames(h3k27ac_umap)),], method = "Euclidean")
dis_blastocyst <- proxy::dist(h3k27ac_umap[grep("blastocyst|E3.5|E3.0", rownames(h3k27ac_umap)),], method = "Euclidean")

max <- median(dis_zygote)
dat <- data.frame(group=c(rep("zygote", length(dis_zygote)), rep("2cell", length(dis_2cell)), rep("4cell", length(dis_4cell)), rep("8cell", length(dis_8cell)), rep("morula", length(dis_morula)), rep("blastocyst", length(dis_blastocyst))),
                  distance=c(dis_zygote/max, dis_2cell/max, dis_4cell/max, dis_8cell/max, dis_morula/max, dis_blastocyst/max))
dat$group <- factor(dat$group, levels = c("zygote", "2cell", "4cell", "8cell", "morula", "blastocyst"))

median <- data.frame(median=c(median(dis_zygote/max), median(dis_2cell/max), median(dis_4cell/max), median(dis_8cell/max), median(dis_morula/max), median(dis_blastocyst/max)),
                     group=c("zygote", "2cell", "4cell", "8cell", "morula", "blastocyst"))

pdf("h3k27ac_eu.pdf")
ggplot(dat,mapping = aes(x = group, y = distance))+
  geom_boxplot(aes(fill = group),width = 0.6)+
  theme_bw()+
  geom_text(data=median, aes(x=group, y=median, label=median), nudge_y=1)
dev.off()

