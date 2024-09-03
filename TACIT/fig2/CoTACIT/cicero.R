setwd("./cicero_h3k4me3_h3k27ac/")

## pseudotime analysis
library(cicero)
library(monocle3)
library(SeuratWrappers)
library(reshape2)
library(cisTopic)
library(Seurat)
library(Signac)
library(tidyr)
set.seed(225)


pathToBams1 <- c('./merge_h3k4me3_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("22114",bamFiles1)]
regions <- "./h3k4me3_h3k27ac_distal_merge.sorted.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac')
c <- as.matrix(cisTopicObject@count.matrix) 

merge <- CreateSeuratObject(counts = c, assay = 'peaks', project = 'merge', min.cells = 1)
merge <- RunTFIDF(merge)
merge <- FindTopFeatures(merge, min.cutoff = 'q0')
merge <- RunSVD(object = merge, assay = 'peaks', reduction.key = 'LSI_', reduction.name = 'lsi')
pdf("./early2cell_merge.pdf")
DepthCor(merge)
merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15),n.neighbors = 12L)
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.3,verbose = FALSE)
DimPlot(object = merge, label = TRUE, pt.size = 2) + NoLegend()
dev.off()

coTarget <- readRDS("../early2cell_coTarget_2.15.rds")
merge@reductions[["umap"]]@cell.embeddings <- coTarget@reductions[["umap.h3k27ac"]]@cell.embeddings
colnames(x = merge[["umap"]]@cell.embeddings) <- paste0("UMAP_", 1:2)
pdf("k27ac.umap.pdf")
DimPlot(object = merge, label = TRUE, pt.size = 2) 
dev.off()
umap <- as.data.frame(merge@reductions[["umap"]]@cell.embeddings)
umap[which(umap$UMAP_1 <0),3] <- "C1"
umap[which(umap$UMAP_1 >0),3] <- "C2"
merge@meta.data$C1C2.ident <- as.factor(umap$V3)
C1 <- subset(merge,cells = colnames(merge[,which(merge@meta.data$C1C2.ident=="C1")]))
C2 <- subset(merge,cells = colnames(merge[,which(merge@meta.data$C1C2.ident=="C2")]))


C1.cds <- as.cell_data_set(C1)
C1.cds <- cluster_cells(cds = C1.cds, reduction_method = "UMAP")
C1.cds <- learn_graph(C1.cds, use_partition = F)

C2.cds <- as.cell_data_set(C2)
C2.cds <- cluster_cells(cds = C2.cds, reduction_method = "UMAP")
C2.cds <- learn_graph(C2.cds, use_partition = F)
# order cells
C1.cds@int_colData$reducedDims$UMAP <- Embeddings(C1, reduction = "umap")
C2.cds@int_colData$reducedDims$UMAP <- Embeddings(C2, reduction = "umap")
C1.cds <- order_cells(C1.cds)
C2.cds <- order_cells(C2.cds)
pdf("pseudotime.pdf")
plot_cells(C1.cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE, 
           graph_label_size=3,
           cell_size = 1.2)
plot_cells(C2.cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE, 
           graph_label_size=3,
           cell_size = 1.2)
dev.off()

# co-accessible
C1_2.cds <- as.CellDataSet(x = C1)
C1.cicero <- make_cicero_cds(C1_2.cds, k=20,reduced_coordinates = C1@reductions[["umap"]]@cell.embeddings)
C2_2.cds <- as.CellDataSet(x = C2)
C2.cicero <- make_cicero_cds(C2_2.cds, k=20,reduced_coordinates = C2@reductions[["umap"]]@cell.embeddings)

# run cicero
C1_conns <- run_cicero(C1.cicero, genomic_coords = './mm10.chrom.sizes', sample_num = 100)
head(C1_conns)
C1_ccans <- generate_ccans(C1_conns)
write.csv(C1_conns,"C1_conns.csv")

C2_conns <- run_cicero(C2.cicero, genomic_coords = './mm10.chrom.sizes', sample_num = 100)
head(C2_conns)
C2_ccans <- generate_ccans(C2_conns)
write.csv(C2_conns,"C2_conns.csv")


######################## choose ZGA-related pairs ######################################
ZGA <- read.csv("./ZGA_up_genes.bed", sep = "\t", header = F)
ZGA$up <- ZGA$V2-2000
ZGA$down <- ZGA$V2+2000
write.table(ZGA[,c(1,5,6)], "ZGA_tss_2kb.bed", sep = "\t", quote = F, col.names = F, row.names = F)

peaks <- as.data.frame(unique(C1_conns$Peak1))
colnames(peaks) <- "peaks"
peaks <- separate(peaks, peaks, c("min", "hh"), sep = ":")
peaks <- separate(peaks, hh, c("hh", "hhh"), sep = "-")
write.table(peaks[,c(1:3)], "peaks.bed", sep = "\t", quote = F, col.names = F, row.names = F)
#bedtools intersect -a peaks.bed -b ZGA_tss_2kb.bed -wa >ZGA_tss_2kb_peaks.bed


C1_conns <- read.csv("C1_conns.csv", sep = ",")
C2_conns <- read.csv("C2_conns.csv", sep = ",")
ZGA_inte <- read.csv("./ZGA_tss_2kb_peaks.bed", sep = "\t", header = F)
ZGA_inte <- unite(ZGA_inte, min, c("V1", "V2"), sep = ":")
ZGA_inte <- unite(ZGA_inte, loci, c("min", "V3"), sep = "-")

C1_conns_int <- C1_conns[which(C1_conns$Peak1 %in% ZGA_inte$loci) ,]
C2_conns_int <- C2_conns[which(C2_conns$Peak1 %in% ZGA_inte$loci) ,]
C1_conns_int$C2_coons <- C2_conns_int$coaccess
C1_conns_int$FC <- (C1_conns_int$coaccess)/(C1_conns_int$C2_coons +0.1)
conns_low <- C1_conns_int[which(C1_conns_int$coaccess < 0 & C1_conns_int$C2_coons < 0),]
conns <- C1_conns_int[!(rownames(C1_conns_int) %in% rownames(conns_low)),]
conns <- na.omit(conns)
conns <- conns[which(conns$coaccess >= 0.2 | conns$C2_coons >=0.2),]
write.table(conns, "C1C2_conns_cutoff0.1.txt", sep = "\t", quote = F)

# Download the GTF associated with this data (mm9) from ensembl and load it
library(rtracklayer)
# download and unzip
# temp <- tempfile()
# download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz", temp)
gene_anno <- readGFF("/media/helab/data1/min/00_reference/Mus_musculus.GRCm38.95.gtf")
# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

pdf("Dppa2_Dppa4_conns_view.pdf")
plot_connections(C1_conns, "chr16", 48119348, 48559442, 
                 gene_model = gene_anno, 
                 viewpoint = c("chr16:48297947-48301007"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
plot_connections(C2_conns, "chr16", 48119348, 48559442, 
                 gene_model = gene_anno, 
                 viewpoint = c("chr16:48297947-48301007"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
dev.off()

pdf("Zscan_conns.pdf")
plot_connections(C1_conns, "chr7", 10495543, 11587976, 
                 gene_model = gene_anno, 
                 #viewpoint = c("chr16:48297947-48301007"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
plot_connections(C2_conns, "chr7", 10495543, 11587976, 
                 gene_model = gene_anno, 
                 #viewpoint = c("chr16:48297947-48301007"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
dev.off()