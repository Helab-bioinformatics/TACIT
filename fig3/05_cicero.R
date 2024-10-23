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
merge <- RunUMAP(object = merge,reduction = 'lsi',dims = c(2:15),n.neighbors = 12L)#或可调试reduction
merge <- FindNeighbors(object = merge,reduction = 'lsi',dims = c(2:15))
merge <- FindClusters(object = merge,algorithm = 3,resolution = 0.3,verbose = FALSE)#或可调试reduction和algorithm
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
gene_anno <- readGFF("./Mus_musculus.GRCm38.95.gtf")
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

####################### promoter and enhancer pairs ##############################
setwd("./cicero_h3k4me3_h3k27ac/ZGA_related/")
conns <- read.csv("C1C2_conns_cutoff0.2.txt", sep = "\t")
h3k4me3 <- read.csv("../h3k4me3_early2cell-q0.05.5kbMerge.bed", sep = "\t", header = F)
h3k4me3$marker <- "h3k4me3"
h3k4me3 <- unite(h3k4me3, min, c("V1", "V2"), sep = ":", remove = F)
h3k4me3 <- unite(h3k4me3, loci, c("min", "V3"), sep = "-", remove = F)
h3k27ac <- read.csv("../h3k27ac_early2cell.distal.peak.bed", sep = "\t", header = F)
h3k27ac$marker <- "h3k27ac"
h3k27ac <- unite(h3k27ac, min, c("V1", "V2"), sep = ":", remove = F)
h3k27ac <- unite(h3k27ac, loci, c("min", "V3"), sep = "-", remove = F)

C1_conns_int <- conns
C1_conns_int[which(C1_conns_int$Peak1 %in% h3k4me3$loci),7] <- 0#0 means promoter, 1 means enhancer
C1_conns_int[which(C1_conns_int$Peak1 %in% h3k27ac$loci),7] <- 1
C1_conns_int[which(C1_conns_int$Peak2 %in% h3k4me3$loci),8] <- 0
C1_conns_int[which(C1_conns_int$Peak2 %in% h3k27ac$loci),8] <- 1
colnames(C1_conns_int) <- c(colnames(C1_conns_int[,1:6]), "Peak1_marker", "Peak2_marker")
C1_conns_int <- C1_conns_int[which(C1_conns_int$Peak1_marker != C1_conns_int$Peak2_marker),]

#################### for promoter-enhancer pairs #####################
#conns <- C1_conns_int[which(C1_conns_int$coaccess >= 0.1 | C1_conns_int$C2_coons >=0.1),]
promoters <- C1_conns_int[which(C1_conns_int$Peak1_marker == 0),]
enhancers <- C1_conns_int[which(C1_conns_int$Peak1_marker == 1),]
colnames(enhancers) <- c("X", "Peak2", "Peak1", "coaccess", "C2_coons", "FC", "Peak2_marker", "Peak1_marker")

promoters <- rbind(promoters, enhancers)
promoters <- separate(promoters, Peak1, c("Peak1_h", "Peak1_hh"), sep = ":", remove = F)
promoters <- separate(promoters, Peak1_hh, c("Peak1_hh", "Peak1_hhh"), sep = "-", remove = F)
promoters <- separate(promoters, Peak2, c("Peak2_h", "Peak2_hh"), sep = ":", remove = F)
promoters <- separate(promoters, Peak2_hh, c("Peak2_hh", "Peak2_hhh"), sep = "-", remove = F)
write.table(promoters[,c("Peak1_h", "Peak1_hh", "Peak1_hhh")], "ZGA_inte_conns_peak1.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(promoters[,c("Peak2_h", "Peak2_hh", "Peak2_hhh")], "ZGA_inte_conns_peak2.bed", sep = "\t", quote = F, col.names = F, row.names = F)

pathToBams1 <- c('../../01_h3k4me3/rename/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("22114",bamFiles1)]
regions <- "ZGA_inte_conns_peak1.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac', min.cells = 0, min.regions = 0)
peak1 <- as.matrix(cisTopicObject@count.matrix) 
pathToBams1 <- c('../../02_h3k27ac/rename/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("22114",bamFiles1)]
regions <- "ZGA_inte_conns_peak2.bed"
cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac',min.cells = 0, min.regions = 0)
peak2 <- as.matrix(cisTopicObject@count.matrix) 

##
peak1_mat <- matrix(1:310770, ncol=90)#3453
peak1_mat[,1] <- promoters$Peak1
for (i in 1:nrow(peak1_mat)) {
  peak1_mat[i,2:90] <- peak1[match(peak1_mat[i,1], rownames(peak1)),]
}
peak1_mat <- as.data.frame(peak1_mat)
colnames(peak1_mat) <- c("loci", colnames(peak1))
peak1_mat$C1conn <- promoters$coaccess
peak1_mat$C2conn <- promoters$C2_coons
peak1_mat$Peak1 <- promoters$Peak1
peak1_mat$Peak2 <- promoters$Peak2
peak1_mat$FC <- abs(peak1_mat$C2conn)/abs(peak1_mat$C1conn)+0.1
up <- peak1_mat[which(peak1_mat$C2conn >peak1_mat$C1conn & peak1_mat$FC >2),]
down <- peak1_mat[which(peak1_mat$C2conn <peak1_mat$C1conn & peak1_mat$FC <0.5),]
noChange <- peak1_mat[!(rownames(peak1_mat) %in% rownames(up)),]
noChange <- noChange[!(rownames(noChange) %in% rownames(down)),]

pseudo_order <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/03_stage_umap/early2cell/h3k27ac_ID.txt", sep = " ")
pseudo_order <- pseudo_order[grep("22114", rownames(pseudo_order)),]
pseudo_order <- separate(pseudo_order, ID, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"), sep = "_")
pseudo_order <- unite(pseudo_order, ID2, c("V1", "V2", "V3", "V4"), sep = "_")
pseudo_order$ID <- paste0(pseudo_order$ID2, ".bam")
pseudo_order <- pseudo_order[order(pseudo_order$C1C2, pseudo_order$pseudo.time),]

up <- up[,c("loci",pseudo_order$ID,colnames(up[,90:95]))]
down <- down[,c("loci",pseudo_order$ID, colnames(down[,90:95]))]
noChange <- noChange[,c("loci",pseudo_order$ID, colnames(noChange[,90:95]))]

##################### for up enhancer-promoter pairs ####################
######## h3k4me3 ##########
up_reads <- up[,2:90]
up_reads <- apply(up_reads, 2, as.numeric)
rownames(up_reads) <- rownames(up)
dat3 <- matrix(1:105087,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat3[,n] <- rowSums(up_reads[,i:m])
  n=n+1
}
rownames(dat3) <- rownames(up)
library(pheatmap)
pdf("up_promoter_enhancer_h3k4me3_3.pdf")
p <- pheatmap(dat3,  
         cluster_rows =T, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D2",
         cellheight=0.2,
         color = colorRampPalette(c("#3c6aa8", "#377bbc","#27b4a3","#c7bd5e","#efdf31"))(500))
p
dev.off()

order <- rownames(up)[p[["tree_row"]][["order"]]]
up <- up[order,]
write.table(up, "up_promoter_enhancer_pairs_h3k4me3_ordered.txt", sep = "\t", quote = F)
# order <- up$loci[p[["tree_row"]][["order"]]]
# write.table(order, "up_promoter_enhancer_pairs_order.txt", quote = F, col.names = F, row.names = F)

############## h3k27ac ####################
peak2_mat <- matrix(1:310770, ncol=90)
peak2_mat[,1] <- as.character(promoters$Peak2)
for (i in 1:nrow(peak2_mat)) {
  peak2_mat[i,2:90] <- peak2[match(peak2_mat[i,1], rownames(peak2)),]
}
peak2_mat <- as.data.frame(peak2_mat)
colnames(peak2_mat) <- c("loci", colnames(peak2))
peak2_mat <- peak2_mat[order,]
peak2_mat <- peak2_mat[,c("loci", pseudo_order$ID)]
up_reads <- peak2_mat[,2:90]
up_reads <- apply(up_reads, 2, as.numeric)
dat3 <- matrix(1:105087,ncol=69)
n=1
for (i in 1:69) {
  m=i+20
  dat3[,n] <- rowSums(up_reads[,i:m])
  n=n+1
}
dat3 <- t(dat3)
dat3 <- scale(dat3)
dat3 <- t(dat3)
dat3[dat3>4] <- 4
pdf("up_promoter_enhancer_h3k27ac_3.pdf")
pheatmap(dat3, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()
write.table(peak2_mat, "up_promoter_enhancer_pairs_h3k27ac.txt", sep = "\t", quote = F)

########## co-binding score ############
score <- up[,c("C1conn", "C2conn")]
score <- score[order,]
pdf("up_promoter_enhancer_cicero.pdf")
pheatmap(score, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#58539f", "#bbbbd6", "#eebabb", "#d86967"))(500))
dev.off()



##################### for down enhancer-promoter pairs ####################
down <- read.csv("down_promoter_enhancer_pairs_h3k4me3.txt", sep = "\t")
down <- down[,c("loci",pseudo_order$ID, colnames(down[,90:95]))]
######## h3k4me3 ##########
down_reads <- down[,2:90]
down_reads <- apply(down_reads, 2, as.numeric)
rownames(down_reads) <- rownames(down)
dat3 <- matrix(1:60375,ncol=69)#875
n=1
for (i in 1:69) {
  m=i+20
  dat3[,n] <- rowSums(down_reads[,i:m])
  n=n+1
}
rownames(dat3) <- rownames(down)
library(pheatmap)
pdf("down_promoter_enhancer_h3k4me3_3.pdf")
p <- pheatmap(dat3, 
              cluster_rows =T, cluster_cols = F, 
              show_rownames = F, scale ="row",
              show_colnames = T, 
              clustering_distance_rows = "canberra",
              clustering_method = "ward.D2",
              cellheight=0.2,
              color = colorRampPalette(c("#3c6aa8", "#377bbc","#27b4a3","#c7bd5e","#efdf31"))(500))
p
dev.off()

order <- rownames(down)[p[["tree_row"]][["order"]]]
down <- down[order,]
write.table(down, "down_promoter_enhancer_pairs_h3k4me3_ordered.txt", sep = "\t", quote = F)

# order <- down$loci[p[["tree_row"]][["order"]]]
# write.table(order, "down_promoter_enhancer_pairs_order.txt", quote = F, col.names = F, row.names = F)

############## h3k27ac ####################
peak2_mat <- matrix(1:310770, ncol=90)
peak2_mat[,1] <- as.character(promoters$Peak2)
for (i in 1:nrow(peak2_mat)) {
  peak2_mat[i,2:90] <- peak2[match(peak2_mat[i,1], rownames(peak2)),]
}
peak2_mat <- as.data.frame(peak2_mat)
colnames(peak2_mat) <- c("loci", colnames(peak2))
####
peak2_mat <- peak2_mat[order,]
peak2_mat <- peak2_mat[,c("loci", pseudo_order$ID)]
down_reads <- peak2_mat[,2:90]
down_reads <- apply(down_reads, 2, as.numeric)
dat3 <- matrix(1:60375,ncol=69)
n=1
for (i in 1:69) {
  m=i+20
  dat3[,n] <- rowSums(down_reads[,i:m])
  n=n+1
}
dat3 <- t(dat3)
dat3 <- scale(dat3)
dat3 <- t(dat3)
dat3[dat3>4] <- 4
pdf("down_promoter_enhancer_h3k27ac_3.pdf")
pheatmap(dat3, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()
write.table(peak2_mat, "down_promoter_enhancer_pairs_h3k27ac.txt", sep = "\t")


########## co-binding score ############
score <- down[,c("C1conn", "C2conn")]
score <- score[order,]
pdf("down_promoter_enhancer_cicero.pdf")
pheatmap(score, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#58539f", "#bbbbd6", "#eebabb", "#d86967"))(500))
dev.off()


##################### for noChange enhancer-promoter pairs ####################
######## h3k4me3 ##########
noChange <- read.csv("noChange_promoter_enhancer_pairs_h3k4me3.txt", sep = "\t")
noChange <- noChange[,c("loci",pseudo_order$ID, colnames(noChange[,90:95]))]

noChange_reads <- noChange[,2:90]
noChange_reads <- apply(noChange_reads, 2, as.numeric)
rownames(noChange_reads) <- rownames(noChange)
dat3 <- matrix(1:72795,ncol=69)#1055
n=1
for (i in 1:69) {
  m=i+20
  dat3[,n] <- rowSums(noChange_reads[,i:m])
  n=n+1
}
rownames(dat3) <- rownames(noChange)
library(pheatmap)
pdf("noChange_promoter_enhancer_h3k4me3_3.pdf")
p <- pheatmap(dat3, 
              cluster_rows =T, cluster_cols = F, 
              show_rownames = F, scale ="row",
              show_colnames = T, 
              clustering_distance_rows = "canberra",
              clustering_method = "ward.D2",
              cellheight=0.2,
              color = colorRampPalette(c("#3c6aa8", "#377bbc","#27b4a3","#c7bd5e","#efdf31"))(500))
p
dev.off()

order <- rownames(noChange)[p[["tree_row"]][["order"]]]
noChange <- noChange[order,]
write.table(noChange, "noChange_promoter_enhancer_pairs_h3k4me3_ordered.txt", sep = "\t", quote = F)

# order <-noChange$loci[p[["tree_row"]][["order"]]]
# write.table(order, "noChange_promoter_enhancer_pairs_order.txt", quote = F, col.names = F, row.names = F)

############## h3k27ac ####################
peak2_mat <- matrix(1:310770, ncol=90)
peak2_mat[,1] <- as.character(promoters$Peak2)
for (i in 1:nrow(peak2_mat)) {
  peak2_mat[i,2:90] <- peak2[match(peak2_mat[i,1], rownames(peak2)),]
}
peak2_mat <- as.data.frame(peak2_mat)
colnames(peak2_mat) <- c("loci", colnames(peak2))
####
peak2_mat <- peak2_mat[order,]
peak2_mat <- peak2_mat[,c("loci", pseudo_order$ID)]
noChange_reads <- peak2_mat[,2:90]
noChange_reads <- apply(noChange_reads, 2, as.numeric)
dat3 <- matrix(1:72795,ncol=69)
n=1
for (i in 1:69) {
  m=i+20
  dat3[,n] <- rowSums(noChange_reads[,i:m])
  n=n+1
}
dat3 <- t(dat3)
dat3 <- scale(dat3)
dat3 <- t(dat3)
dat3[dat3>4] <- 4
pdf("noChange_promoter_enhancer_h3k27ac_3.pdf")
pheatmap(dat3, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()
write.table(peak2_mat, "noChange_promoter_enhancer_pairs_h3k27ac.txt", sep = "\t")

########## co-binding score ############
score <- noChange[,c("C1conn", "C2conn")]
score <- score[order,]
score[nrow(score)+1,] <- c(0.85, -0.4)
pdf("noChange_promoter_enhancer_cicero.pdf")
pheatmap(score, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#58539f", "#bbbbd6", "#eebabb", "#d86967"))(500))
dev.off()


##################################### RNA ##################################################

############################## integrate h3k4me3 with RNA ##############################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/coTarget_early2cell")
h3k4me3 <- readRDS("h3k4me3_cotarget_early2cell.rds")
h3k4me3_counts <- as.data.frame(h3k4me3@assays[["peaks"]]@counts)
geneActivity <- read.csv("./h3k4me3_early2cell_5kb_Activityrix.txt", sep = " ",check.names=F)
h3k4me3[["ACTIVITY"]] <- CreateAssayObject(counts = geneActivity)
DefaultAssay(h3k4me3) <- "ACTIVITY"
h3k4me3 <- FindVariableFeatures(h3k4me3, nfeatures = 10000)
h3k4me3 <- NormalizeData(h3k4me3)
h3k4me3 <- ScaleData(h3k4me3)

RNA <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/07_repeat/20221107_q10/integration/ZGA_rna.rds")
RNA <- subset(x=RNA, idents =c("early2cell", "mid2cell", "late2cell", "zygote"))

transfer.anchors <- FindTransferAnchors(reference = RNA,
                                        query = h3k4me3,
                                        features = VariableFeatures(object = RNA),
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)
#Annotate scATAC-seq cells via label transfer
cellID.predictions <- TransferData(anchorset = transfer.anchors, refdata = RNA$ID,
                                   weight.reduction = h3k4me3[["lsi"]], dims = 2:ncol(h3k4me3[["lsi"]]), k.weight = 50)
cellID.predictions = subset(cellID.predictions, select = c('predicted.id', 'prediction.score.max'))
names(cellID.predictions)[1] = 'Predicted_Cell_ID'
h3k4me3 <- AddMetaData(h3k4me3, metadata = cellID.predictions)
h3k4me3_ID <- as.data.frame(h3k4me3@meta.data[["Predicted_Cell_ID"]])

##############################################################################################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/coTarget_early2cell/cicero")
up <- read.csv("up_promoter_enhancer_pairs_order.txt", sep = "\t", header = F)
down <- read.csv("down_promoter_enhancer_pairs_order.txt", sep = "\t", header = F)
noChange <- read.csv("noChange_promoter_enhancer_pairs_order.txt", sep = "\t", header = F)

library(tidyr)
noChange <- separate(noChange, V1, c("min", "hh"), sep = ":")
noChange <- separate(noChange, hh, c("hh", "hhh"), sep = "-")
write.table(noChange, "noChange_promoter_enhancer_pairs_order.bed", sep = "\t", quote = F, col.names = F, row.names = F)

counts <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/GSE45719_scRNA.csv")
counts <- counts[!duplicated(counts$Gene_symbol),]
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[,grep("zygote|2cell", colnames(counts))]

pseudo_order <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/03_stage_umap/early2cell/h3k27ac_ID.txt", sep = " ")
pseudo <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/03_stage_umap/early2cell/h3k27ac_pseudotime.txt", sep = "\t")
pseudo_order$chip.pseuotime <- pseudo$pseudo_order[match(rownames(pseudo_order), rownames(pseudo))]
pseudo_order <- pseudo_order[grep("22114", rownames(pseudo_order)),]
pseudo_order <- separate(pseudo_order, ID, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"), sep = "_")
pseudo_order <- unite(pseudo_order, ID2, c("V1", "V2", "V3", "V4"), sep = "_")
pseudo_order$ID <- paste0(pseudo_order$ID2, ".bam")
pseudo_order <- pseudo_order[order(pseudo_order$C1C2, pseudo_order$chip.pseuotime),]

up_anno <- read.csv("up_promoter_enhancer_pairs_order_anno.txt", sep = "\t")
colnames(up_anno) <- c("order", colnames(up_anno[,2:ncol(up_anno)]))
up_anno <- up_anno[order(up_anno$order),]
RNA <- counts[up_anno$Gene.Name,]
pseudo_order_RNA <- pseudo_order[order(pseudo_order$pseudo.time),]
RNA <- RNA[,pseudo_order_RNA$predict_ID]
pdf("up_promoter_enhancer_RNA_2.pdf")
pheatmap(RNA, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
dev.off()


down_anno <- read.csv("down_promoter_enhancer_pairs_order_anno.txt", sep = "\t")
colnames(down_anno) <- c("order", colnames(down_anno[,2:ncol(down_anno)]))
down_anno <- down_anno[order(down_anno$order),]
RNA <- counts[down_anno$Gene.Name,]
RNA <- RNA[,pseudo_order_RNA$predict_ID]
pdf("down_promoter_enhancer_RNA_2.pdf")
pheatmap(RNA, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
dev.off()

noChange_anno <- read.csv("noChange_promoter_enhancer_pairs_order_anno.txt", sep = "\t")
colnames(noChange_anno) <- c("order", colnames(noChange_anno[,2:ncol(noChange_anno)]))
noChange_anno <- noChange_anno[order(noChange_anno$order),]
RNA <- counts[noChange_anno$Gene.Name,]
RNA <- RNA[,pseudo_order_RNA$predict_ID]
pdf("noChange_promoter_enhancer_RNA_2.pdf")
pheatmap(RNA, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
dev.off()


################################ for co-clustering ################################
setwd("./ZGA_related/")

pseudo_order <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/03_stage_umap/early2cell/h3k27ac_ID.txt", sep = " ")
pseudo_order <- pseudo_order[grep("22114", rownames(pseudo_order)),]
pseudo_order <- separate(pseudo_order, ID, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"), sep = "_")
pseudo_order <- unite(pseudo_order, ID2, c("V1", "V2", "V3", "V4"), sep = "_")
pseudo_order$ID <- paste0(pseudo_order$ID2, ".bam")
pseudo_order <- pseudo_order[order(pseudo_order$C1C2, pseudo_order$pseudo.time),]

################################### up #################################################
up_k4 <- read.csv("./up_promoter_enhancer_pairs_h3k4me3_ordered.txt", sep = "\t")
up_k4 <- up_k4[,c("loci",pseudo_order$ID,colnames(up_k4[,92:96]))]
up_reads <- up_k4[,2:90]
up_reads <- apply(up_reads, 2, as.numeric)
rownames(up_reads) <- rownames(up_k4)
dat1 <- matrix(1:105087,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat1[,n] <- rowSums(up_reads[,i:m])
  n=n+1
}
rownames(dat1) <- rownames(up_k4)

up_k27 <- read.csv("./up_promoter_enhancer_pairs_h3k27ac.txt", sep = "\t")
up_k27 <- up_k27[,c("loci",pseudo_order$ID)]
up_reads <- up_k27[,2:90]
up_reads <- apply(up_reads, 2, as.numeric)
rownames(up_reads) <- rownames(up_k27)
dat2 <- matrix(1:105087,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat2[,n] <- rowSums(up_reads[,i:m])
  n=n+1
}
rownames(dat2) <- rownames(up_k27)

dat1 <- t(dat1)
dat1 <- scale(dat1)
dat1 <- t(dat1)
dat1[dat1>4] <- 4
dat2 <- t(dat2)
dat2 <- scale(dat2)
dat2 <- t(dat2)
dat2[dat2>4] <- 4

dat3 <- cbind(dat1, dat2)
pdf("up_promoter_enhancer_h3k4me3_h3k27ac_cocluster.pdf")
p <- pheatmap(dat3, 
              cluster_rows =T, cluster_cols = F, 
              show_rownames = F, scale ="none",
              show_colnames = T, 
              cellheight=0.2, clustering_method = "complete", 
              clustering_distance_rows = "euclidean",
              color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
p
dev.off()

order <- up_k4$loci[p[["tree_row"]][["order"]]]
order <- as.data.frame(order)
order <- separate(order, order, c("min", "h"), sep = ":")
order <- separate(order, h, c("h", "hh"), sep = "-")
write.table(order, "up_coclustering_ordered.txt", sep = "\t", quote = F, col.names = F, row.names = F)

################## score ###################
row <- rownames(up_k4)[p[["tree_row"]][["order"]]]
up_k4 <- up_k4[row,]
score <- up_k4[,c("C1conn", "C2conn")]
score[nrow(score)+1,] <- c(0.9, -0.9)
pdf("up_promoter_enhancer_cicero.pdf")
pheatmap(score, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#58539f", "#bbbbd6", "#eebabb", "#d86967"))(500))
dev.off()

dat2 <- dat2[row,]
pdf("up_promoter_enhancer_h3k27ac.pdf")
pheatmap(dat2, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()

################### RNA ####################
counts <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/GSE45719_scRNA.csv")
counts <- counts[!duplicated(counts$Gene_symbol),]
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[,grep("zygote|2cell", colnames(counts))]

up_anno <- read.csv("up_coclustering_ordered_anno.txt", sep = "\t")
colnames(up_anno) <- c("order", colnames(up_anno[,2:ncol(up_anno)]))
RNA <- counts[up_anno$Gene.Name,]
pseudo_order_RNA <- pseudo_order[order(pseudo_order$pseudo.time),]
RNA <- RNA[,pseudo_order_RNA$predict_ID]
pdf("up_promoter_enhancer_cocluster_RNA.pdf")
pheatmap(RNA, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
dev.off()


################################## down #####################################
down_k4 <- read.csv("./down_promoter_enhancer_pairs_h3k4me3_ordered.txt", sep = "\t")
down_k4 <- down_k4[,c("loci",pseudo_order$ID,colnames(down_k4[,92:96]))]
down_reads <- down_k4[,2:90]
down_reads <- apply(down_reads, 2, as.numeric)
rownames(down_reads) <- rownames(down_k4)
dat1 <- matrix(1:60375,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat1[,n] <- rowSums(down_reads[,i:m])
  n=n+1
}
rownames(dat1) <- rownames(down_k4)

down_k27 <- read.csv("./down_promoter_enhancer_pairs_h3k27ac.txt", sep = "\t")
down_k27 <- down_k27[,c("loci",pseudo_order$ID)]
down_reads <- down_k27[,2:90]
down_reads <- apply(down_reads, 2, as.numeric)
rownames(down_reads) <- rownames(down_k27)
dat2 <- matrix(1:60375,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat2[,n] <- rowSums(down_reads[,i:m])
  n=n+1
}
rownames(dat2) <- rownames(down_k27)

dat1 <- t(dat1)
dat1 <- scale(dat1)
dat1 <- t(dat1)
dat1[dat1>4] <- 4
dat2 <- t(dat2)
dat2 <- scale(dat2)
dat2 <- t(dat2)
dat2[dat2>4] <- 4

dat3 <- cbind(dat1, dat2)
pdf("down_promoter_enhancer_h3k4me3_h3k27ac_cocluster.pdf")
p <- pheatmap(dat3, 
              cluster_rows =T, cluster_cols = F, 
              show_rownames = F, scale ="none",
              show_colnames = T, 
              cellheight=0.2, clustering_method = "average", 
              clustering_distance_rows = "manhattan",
              color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
p
dev.off()

order <- down_k4$loci[p[["tree_row"]][["order"]]]
order <- as.data.frame(order)
order <- separate(order, order, c("min", "h"), sep = ":")
order <- separate(order, h, c("h", "hh"), sep = "-")
write.table(order, "down_coclustering_ordered.txt", sep = "\t", quote = F, col.names = F, row.names = F)

row <- rownames(down_k4)[p[["tree_row"]][["order"]]]
down_k4 <- down_k4[row,]
score <- down_k4[,c("C1conn", "C2conn")]
score[nrow(score)+1,] <- c(0.9, -0.9)
pdf("down_promoter_enhancer_cicero.pdf")
pheatmap(score, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#58539f", "#bbbbd6", "#eebabb", "#d86967"))(500))
dev.off()

dat2 <- dat2[row,]
pdf("down_promoter_enhancer_h3k27ac.pdf")
pheatmap(dat2, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()

################### RNA ####################
counts <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/GSE45719_scRNA.csv")
counts <- counts[!duplicated(counts$Gene_symbol),]
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[,grep("zygote|2cell", colnames(counts))]

down_anno <- read.csv("down_coclustering_ordered_anno.txt", sep = "\t")
colnames(down_anno) <- c("order", colnames(down_anno[,2:ncol(down_anno)]))
RNA <- counts[down_anno$Gene.Name,]
pseudo_order_RNA <- pseudo_order[order(pseudo_order$pseudo.time),]
RNA <- RNA[,pseudo_order_RNA$predict_ID]
pdf("down_promoter_enhancer_cocluster_RNA.pdf")
pheatmap(RNA, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
dev.off()


################################## noChange #####################################
noChange_k4 <- read.csv("./noChange_promoter_enhancer_pairs_h3k4me3_ordered.txt", sep = "\t")
noChange_k4 <- noChange_k4[,c("loci",pseudo_order$ID,colnames(noChange_k4[,92:96]))]
noChange_reads <- noChange_k4[,2:90]
noChange_reads <- apply(noChange_reads, 2, as.numeric)
rownames(noChange_reads) <- rownames(noChange_k4)
dat1 <- matrix(1:72795,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat1[,n] <- rowSums(noChange_reads[,i:m])
  n=n+1
}
rownames(dat1) <- rownames(noChange_k4)

noChange_k27 <- read.csv("./noChange_promoter_enhancer_pairs_h3k27ac.txt", sep = "\t")
noChange_k27 <- noChange_k27[,c("loci",pseudo_order$ID)]
noChange_reads <- noChange_k27[,2:90]
noChange_reads <- apply(noChange_reads, 2, as.numeric)
rownames(noChange_reads) <- rownames(noChange_k27)
dat2 <- matrix(1:72795,ncol=69)#1523
n=1
for (i in 1:69) {
  m=i+20
  dat2[,n] <- rowSums(noChange_reads[,i:m])
  n=n+1
}
rownames(dat2) <- rownames(noChange_k27)

dat1 <- t(dat1)
dat1 <- scale(dat1)
dat1 <- t(dat1)
dat1[dat1>4] <- 4
dat2 <- t(dat2)
dat2 <- scale(dat2)
dat2 <- t(dat2)
dat2[dat2>4] <- 4

dat3 <- cbind(dat1, dat2)
pdf("noChange_promoter_enhancer_h3k4me3_h3k27ac_cocluster.pdf")
p <- pheatmap(dat3, 
              cluster_rows =T, cluster_cols = F, 
              show_rownames = F, scale ="none",
              show_colnames = T, 
              cellheight=0.2,
              clustering_distance_rows="euclidean",	
              clustering_method="ward.D",
              color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
p
dev.off()

order <- noChange_k4$loci[p[["tree_row"]][["order"]]]
order <- as.data.frame(order)
order <- separate(order, order, c("min", "h"), sep = ":")
order <- separate(order, h, c("h", "hh"), sep = "-")
write.table(order, "noChange_coclustering_ordered.txt", sep = "\t", quote = F, col.names = F, row.names = F)

row <- rownames(noChange_k4)[p[["tree_row"]][["order"]]]
noChange_k4 <- noChange_k4[row,]
score <- noChange_k4[,c("C1conn", "C2conn")]
score[nrow(score)+1,] <- c(0.9, -0.9)
pdf("noChange_promoter_enhancer_cicero.pdf")
pheatmap(score, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#58539f", "#bbbbd6", "#eebabb", "#d86967"))(500))
dev.off()


dat2 <- dat2[row,]
pdf("noChange_promoter_enhancer_h3k27ac.pdf")
pheatmap(dat2, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()

################### RNA ####################
counts <- read.csv("/media/helab/data1/min/00_reference/mouse_embryo/genes/scRNA-seq_2014/GSE45719_scRNA.csv")
counts <- counts[!duplicated(counts$Gene_symbol),]
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[,grep("zygote|2cell", colnames(counts))]

noChange_anno <- read.csv("noChange_coclustering_ordered_anno.txt", sep = "\t")
colnames(noChange_anno) <- c("order", colnames(noChange_anno[,2:ncol(noChange_anno)]))
RNA <- counts[noChange_anno$Gene.Name,]
pseudo_order_RNA <- pseudo_order[order(pseudo_order$pseudo.time),]
RNA <- RNA[,pseudo_order_RNA$predict_ID]
pdf("noChange_promoter_enhancer_cocluster_RNA.pdf")
pheatmap(RNA, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="row",
         show_colnames = T, 
         cellheight=0.2,
         color = colorRampPalette(c("#3b6bab","#6bb2dc","#e6d28d","#ac2728"))(500))
dev.off()
