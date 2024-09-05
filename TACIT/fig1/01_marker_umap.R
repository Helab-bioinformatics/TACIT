rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999) 
options(stringsAsFactors = F) 

#cisTopic for countMatrix
library(cisTopic)
library(ggplot2)
library(ggpubr)
library(ggsignif)

pathToBams1 <- c('./allBam/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("txt",invert = T, bamFiles1)]

regions <- "./mm10.5K.windows.bed"
#regions <- "h3k4me1.peak-0.05_peaks.broadPeak"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k4me1')
c <- as.matrix(cisTopicObject@count.matrix) 

#signac
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
library(stringr)

h3k4me1 <- CreateSeuratObject(
  counts = c,
  assay = 'peaks',
  project = 'h3k4me1',
  min.cells = 1
)
h3k4me1

pdf("h3k4me1_9.7.pdf")
h3k4me1 <- RunTFIDF(h3k4me1)
h3k4me1 <- FindTopFeatures(h3k4me1, min.cutoff = '5')
h3k4me1 <- RunSVD(
  object = h3k4me1,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)
DepthCor(h3k4me1)

h3k4me1 <- RunUMAP(
  object = h3k4me1,
  reduction = 'lsi',
  dims = c(2:15)
)#或可调试reduction
h3k4me1 <- FindNeighbors(
  object = h3k4me1,
  reduction = 'lsi',
  dims = c(2:15)
)
h3k4me1 <- FindClusters(
  object = h3k4me1,
  algorithm = 1,
  resolution = 0.5,
  verbose = FALSE
)#或可调试reduction和algorithm
DimPlot(object = h3k4me1, label = TRUE, pt.size = 2) + NoLegend()

label <- rownames(as.data.frame(h3k4me1@active.ident))
label <- as.data.frame(label)
label[grep("zygote",label[,1]),2] <- "zygote"
label[grep("2cell",label[,1]),2] <- "2cell"
label[grep("4cell",label[,1]),2] <- "4cell"
label[grep("E2.0|8cell",label[,1]),2] <- "8cell"
label[grep("E2.5|morula",label[,1]),2] <- "morula"
label[grep("E3.0|E3.5|blasto",label[,1]),2] <- "E3.0"

h3k4me1@meta.data$group.ident <- as.factor(label[,2])
pdf("h3k4me1.5kb_noOuter.pdf")
DimPlot(object = h3k4me1, group.by = "group.ident" ,label = TRUE, pt.size = 2) 
DimPlot(object = h3k4me1 ,label = TRUE,cols = DiscretePalette(15, palette = NULL),pt.size = 2) 

#depth
depth <- read.csv("./allBam/Min/depth.txt", sep = "\t", header = F)
# depth <- separate(depth, ID, c("1", "2", "3", "4", "5"), sep = "_")
# depth <- unite(depth, "cell", 1,2,3,4,5, sep=".")
#depth[,5] <- paste0(depth[,1], "_rmdup.bam")
row.names(depth) <- depth[,1]
correctRaw <- label$label
depth <- depth[correctRaw,]
depth[,3] <- log10(depth[,2])

h3k4me1@meta.data$depth.ident <- as.numeric(depth[,2])
h3k4me1@meta.data$log10.depth.ident <- as.numeric(depth[,3])

FeaturePlot(
  object = h3k4me1,
  features = 'depth.ident',
  pt.size = 1.5,
  max.cutoff = 'q95',
  min.cutoff = 'q5'
) 
FeaturePlot(
  object = h3k4me1,
  features = 'log10.depth.ident',
  pt.size = 1.5,
  max.cutoff = 'q95',
  min.cutoff = 'q5'
) 
dev.off()

write.table(label, "label.txt", quote = F)
save(c, h3k4me1, label, file = "./h3k4me1.Rdata")
