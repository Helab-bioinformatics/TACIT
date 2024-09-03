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

# 对Circos.Thickness进行归一化
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


############################## overlap with MERVL-related H3K4me3-H3K27ac interactions ##################################
setwd("./inte_cicero/")
interactions <- read.table("../late2cell_50kb_interactions.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
interactions <- interactions[which(interactions$FDR.Benjamini..based.on.1.84e.09.total.tests. <0.00001),]
interactions <- interactions[which(interactions$Interaction.Reads > 100),]
write.table(interactions[,c("chr.1.", "start.1.", "end.1.", "InteractionID")], "./late2cell_50kb_interactions_reads100_peak1.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(interactions[,c("chr.2.", "start.2.", "end.2.", "InteractionID")], "./late2cell_50kb_interactions_reads100_peak2.bed", sep = "\t", quote = F, row.names = F, col.names = F)

###得到MERVL相关的所有links的情况
MERVL_pr_peak1 <- read.csv("./MERVL_inte_promoter_cutoff0.1_peak1.bed", sep = "\t", header = F)
MERVL_pr_peak2 <- read.csv("./MERVL_inte_promoter_cutoff0.1_peak2.bed", sep = "\t", header = F)
MERVL_pr <- cbind(MERVL_pr_peak1, MERVL_pr_peak2)
MERVL_en_peak1 <- read.csv("./MERVL_inte_enhancer_cutoff0.1_peak1.bed", sep = "\t", header = F)
MERVL_en_peak2 <- read.csv("./MERVL_inte_enhancer_cutoff0.1_peak2.bed", sep = "\t", header = F)
MERVL_en <- cbind(MERVL_en_peak1, MERVL_en_peak2)
MERVL <- rbind(MERVL_pr, MERVL_en)
colnames(MERVL) <- c("chr_1", "start_1", "end_1", "chr_2", "start_2", "end_2")
MERVL$cicerID <- paste0("cicero_", 1:nrow(MERVL))
write.table(MERVL[,c("chr_1", "start_1", "end_1","cicerID")], "./MERVL_alllinks_cutoff0.1_peak1.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(MERVL[,c("chr_2", "start_2", "end_2", "cicerID")], "./MERVL_alllinks_cutoff0.1_peak2.bed", sep = "\t", quote = F, row.names = F, col.names = F)

#######
link1 <- read.csv("./MERVL_alllinks_cutoff0.1_peak1_inte_late2cell_50kb_interactions.bed", sep = "\t", header = F)
link2 <- read.csv("./MERVL_alllinks_cutoff0.1_peak2_inte_late2cell_50kb_interactions.bed", sep = "\t", header = F)
link <- intersect(link1$V4, link2$V4)# 118

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
dat_hic <- data.frame(group=c("HiC-links", "HiC-only"), var=c(4188, 36096))
dat_links <- data.frame(group=c("Links-HiC", "Links-only"), var=c(118, 1070))

library(ggplot2)
library(patchwork)
pdf("link_HiC_ratio.pdf")
p1<- ggplot(dat_hic, aes(x="", y=var, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
p2<- ggplot(dat_links, aes(x="", y=var, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  blank_theme +
  theme(axis.text.x=element_blank())
p1+p2
dev.off()


#################################### 只看ZGA激活后的MERVL links的情况 #############################
promoters <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/coTarget_early2cell/cicero_h3k4me3_h3k27ac/all_tss5kb/promoter_enhancer_pairs_cutoff0.1.txt", sep = "\t")
inte_en <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/coTarget_early2cell/cicero_h3k4me3_h3k27ac/all_tss5kb/C1C2_conns_cutoff0.1_enhancer_inte_MERVL.bed", sep = "\t", header = F)
inte_en <- unite(inte_en, "min", c("V1", "V2"), sep = ":")
inte_en <- unite(inte_en, "loci", c("min", "V3"), sep = "-")
inte_en <- promoters[which(promoters$Peak2 %in% inte_en$loci),]
inte_en <- inte_en[,c("Peak1", "Peak2", "coaccess", "C2_coons", "Peak1_marker", "Peak2_marker")]
inte_en <- separate(inte_en, Peak1, c("Peak1_h", "Peak1_hh"), sep = ":", remove = F)
inte_en <- separate(inte_en, Peak1_hh, c("Peak1_hh", "Peak1_hhh"), sep = "-", remove = F)
inte_en <- separate(inte_en, Peak2, c("Peak2_h", "Peak2_hh"), sep = ":", remove = F)
inte_en <- separate(inte_en, Peak2_hh, c("Peak2_hh", "Peak2_hhh"), sep = "-", remove = F)

inte_pr <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/coTarget_early2cell/cicero_h3k4me3_h3k27ac/all_tss5kb/C1C2_conns_cutoff0.1_promoter_inte_MERVL.bed", sep = "\t", header = F)
inte_pr <- unite(inte_pr, "min", c("V1", "V2"), sep = ":")
inte_pr <- unite(inte_pr, "loci", c("min", "V3"), sep = "-")
inte_pr <- promoters[which(promoters$Peak1 %in% inte_pr$loci),]
inte_pr <- inte_pr[,c("Peak1", "Peak2", "coaccess", "C2_coons", "Peak1_marker", "Peak2_marker")]
inte_pr <- separate(inte_pr, Peak1, c("Peak1_h", "Peak1_hh"), sep = ":", remove = F)
inte_pr <- separate(inte_pr, Peak1_hh, c("Peak1_hh", "Peak1_hhh"), sep = "-", remove = F)
inte_pr <- separate(inte_pr, Peak2, c("Peak2_h", "Peak2_hh"), sep = ":", remove = F)
inte_pr <- separate(inte_pr, Peak2_hh, c("Peak2_hh", "Peak2_hhh"), sep = "-", remove = F)

MERVL_inte <- rbind(inte_en, inte_pr)
MERVL_inte_up <- MERVL_inte[which(MERVL_inte$C2_coons >MERVL_inte$coaccess),]
MERVL_inte_up$cieroID <- paste0("cicero_", 1:nrow(MERVL_inte_up))
write.table(MERVL_inte_up[,c("Peak1_h", "Peak1_hh", "Peak1_hhh")], "MERVL_alllinks_cutoff0.1_up_peak1.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(MERVL_inte_up[,c("Peak2_h", "Peak2_hh", "Peak2_hhh")], "MERVL_alllinks_cutoff0.1_up_peak2.bed", sep = "\t", quote = F, col.names = F, row.names = F)

###########得到的比例更少，故没用######




################################### 统计peak1和peak2都能对应上的 ###################################
mervl_peak1 <- read.csv("./tmp/MERVL_alllinks_cutoff0.1_peak1_50kb.bed", sep = "\t", header = F)
mervl_peak1 <- mervl_peak1[!duplicated(mervl_peak1$V7),]
mervl_peak2 <- read.csv("./tmp/MERVL_alllinks_cutoff0.1_peak2_50kb.bed", sep = "\t", header = F)
mervl_peak2 <- mervl_peak2[!duplicated(mervl_peak2$V7),]
rownames(mervl_peak2) <- mervl_peak2$V7
mervl_peak2 <- mervl_peak2[mervl_peak1$V7,]

mervl_inte <- cbind(mervl_peak1, mervl_peak2)
colnames(mervl_inte) <- c(paste0("peak1_", 1:7), paste0("peak2_", 1:7))
mervl_inte <- unite(mervl_inte, id, c("peak1_1", "peak1_2","peak1_3", "peak2_1", "peak2_2","peak2_3"))

interactions <- read.csv("./late2cell_50kb_interactions_reads100_peak1_inte_MERVL.bed", sep = "\t", header = F)
all_inte <- read.table("../late2cell_50kb_interactions.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
inte <- all_inte[which(all_inte$InteractionID %in% interactions$V4),]
inte <- unite(inte, id, c("chr.1.", "start.1.", "end.1.","chr.2.", "start.2.", "end.2."), sep = "_")

a <- intersect(mervl_inte$id, inte$id)



#######
link1 <- read.csv("./MERVL_alllinks_cutoff0.1_peak1_inte_late2cell_50kb_interactions_reads100.bed", sep = "\t", header = F)
link2 <- read.csv("./MERVL_alllinks_cutoff0.1_peak2_inte_late2cell_50kb_interactions_reads100.bed", sep = "\t", header = F)
link <- intersect(link1$V4, link2$V4)# 118
