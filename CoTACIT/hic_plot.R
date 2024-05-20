setwd("/media/helab/data2/min/HiC_analysis")

library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(GenomicRanges)

interactions <- read.table("./late2cell_50kb_interactions.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

library(tidyr)
inter <- interactions[,c("chr.1.", "start.1.", "end.1.","chr.2.", "start.2.", "end.2.","Circos.Thickness")]
inter <- unite(inter, min, c("chr.1.", "start.1."), sep = ":")
inter <- unite(inter, Peak1, c("min", "end.1."), sep = "-")
inter <- unite(inter, min, c("chr.2.", "start.2."), sep = ":")
inter <- unite(inter, Peak2, c("min", "end.2."), sep = "-")

# normalization
a <- inter$Circos.Thickness
normalized_a <- (a - min(a)) / (max(a) - min(a))
inter[,3] <- normalized_a


library(rtracklayer)
library(cicero)
# download and unzip
# temp <- tempfile()
# download.file("ftp://ftp.ensembl.org/pub/release-65/gtf/mus_musculus/Mus_musculus.NCBIM37.65.gtf.gz", temp)
gene_anno <- readGFF("/media/helab/data1/min/00_reference/Mus_musculus.GRCm38.95.gtf")
# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

# inter <- read.csv("late2cell_interactions_plot.txt", sep = "\t")
# inter_se <- inter[grep("chr10:1172", inter$Peak1),]

colnames(inter) <- c("Peak1","Peak2","coaccess")
pdf("MERVL_cases_view_50kb.pdf")
plot_connections(inter, "chr10", 117228640, 118210310, 
                 gene_model = gene_anno, 
                 viewpoint = c("chr10:117,709,807-117,729,000"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 

plot_connections(inter, "chr11", 61094915, 62011895,
                 gene_model = gene_anno, 
                 viewpoint = c("chr11:61529266-61577543"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
plot_connections(inter, "chr12", 83633216, 84535523,
                 gene_model = gene_anno, 
                 viewpoint = c("chr12:84078212-84090529"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
plot_connections(inter, "chr15", 98300190, 99081677, 
                 gene_model = gene_anno, 
                 viewpoint = c("chr15:98,687,497-98,694,344"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest")
plot_connections(inter, "chr16", 20011038, 20748691,
                 gene_model = gene_anno, 
                 viewpoint = c("chr16:20336660-20423070"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
plot_connections(inter, "chr17", 26934647, 27648110,
                 gene_model = gene_anno, 
                 viewpoint = c("chr17:27286339-27296419"),
                 coaccess_cutoff = 0,
                 connection_width = .5,
                 comparison_connection_width = .5,
                 alpha_by_coaccess = T, 
                 include_axis_track = F,
                 collapseTranscripts = "longest") 
dev.off()
