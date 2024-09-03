setwd("/media/helab/helab-7/LM_all/20240427_IVF_embryo_TACIT/IVF_embryo_TACIT")

####################################################################################################################
######################################### compare ZGA scare between the two cells ##################################
#################### h3k4me3 early #####################
h3k4me3_early <- readRDS("./h3k4me3_early_2kbbin.rds")
label <- rownames(as.data.frame(h3k4me3_early@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k4me3_early@meta.data[["zga.signal.nor"]])

####### 调整ZGA score ######
path <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/test/"
bed <- list.files(path)
k4_early_fc_titration <- vector("list", length(bed))
dif_portion <- c()
for (i in 1:length(bed)) {
  pathToBams1 <- c('./early2cell_h3k4me3/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
  regions <- paste0(path, bed[i])
  cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3', min.cells = 0, min.regions = 0)
  h3k4me3_zga_early <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
  h3k4me3_zga_early <- colSums(h3k4me3_zga_early)
  
  cell_id <- read.csv("./early2cell_k4me3_all_QC.txt", sep = "\t", header = F)
  cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
  h3k4me3_zga_early <- h3k4me3_zga_early[which(names(h3k4me3_zga_early) %in% cell_id$id)]
  h3k4me3_all_counts <- as.numeric(h3k4me3_early@meta.data[["nCount_peaks"]])
  label$zga_signal <- h3k4me3_zga_early/h3k4me3_all_counts
  
  embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
  k4_early_fc <- c()
  for (j in 1:length(embryos)) {
    dat_select <- label[grep(embryos[j], label$x5),]
    k4_early_fc_tmp <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
    k4_early_fc <- c(k4_early_fc, k4_early_fc_tmp)
  }
  k4_early_fc
  k4_early_fc_titration[[i]] <- k4_early_fc
  dif <- k4_early_fc[k4_early_fc > 1.3 | k4_early_fc < 0.8]
  dif_portion_tmp <- length(dif)/length(k4_early_fc)
  dif_portion <- c(dif_portion, dif_portion_tmp)
}
names(dif_portion) <- bed
dif_portion_k4_early <- dif_portion # MajorZGA_genes_tss.bed 或者 majorZGA.genebody.bed 上两个细胞差异最大

###### use MajorZGA_genes_tss.bed #####
label <- rownames(as.data.frame(h3k4me3_early@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k4me3_early@meta.data[["zga.signal.nor"]])

pathToBams1 <- c('./early2cell_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/tss/MajorZGA_genes_tss.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k4me3', min.cells = 0, min.regions = 0)
h3k4me3_zga_early <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
h3k4me3_zga_early <- colSums(h3k4me3_zga_early)

cell_id <- read.csv("./early2cell_k4me3_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k4me3_zga_early <- h3k4me3_zga_early[which(names(h3k4me3_zga_early) %in% cell_id$id)]
h3k4me3_all_counts <- as.numeric(h3k4me3_early@meta.data[["nCount_peaks"]])
label$zga_signal <- h3k4me3_zga_early/h3k4me3_all_counts

embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
dim <- length(embryos) * 4
k4_early_zga <- matrix(1:dim, ncol = 4)
k4_early_zga <- as.data.frame(k4_early_zga)
for (j in 1:length(embryos)) {
  dat_select <- label[grep(embryos[j], label$x5),]
  k4_early_zga[j,1] <- embryos[j]
  k4_early_zga[j,2] <- dat_select$zga_signal[1]
  k4_early_zga[j,3] <- dat_select$zga_signal[2]
  k4_early_zga[j,4] <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
}
colnames(k4_early_zga) <- c("embryo", "cell1", "cell2", "foldChange")
k4_early_zga$low_zga <- pmin(k4_early_zga$cell1, k4_early_zga$cell2)
k4_early_zga$high_zga <- pmax(k4_early_zga$cell1, k4_early_zga$cell2)

library(ggplot2)
# 假设你的数据框是df，有ID列标识每行
df <- data.frame(ID = k4_early_zga$embryo, A = k4_early_zga$low_zga, B = k4_early_zga$high_zga)
# 将数据从宽格式转换为长格式
df_long <- reshape2::melt(df, id.vars = "ID", measure.vars = c("A", "B"))
pdf("./h3k4me3_early_c1_c2_zga_score.pdf")
ggplot(df_long, aes(x = variable, y = value, group = ID)) +
  geom_line(color = "lightgray") +  # 连接线
  geom_point() +  # 可选，增加点
  theme_bw() +  # 选择一个简洁的主题
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1.2) +  # 平均值连线
  stat_compare_means(comparisons = list(c("A", "B"))) +  # 统计显著性
  labs(x = "组别", y = "数值大小", title = "AB两列数值的连接线图")
dev.off()

  
  ########################################### h3k4me3 late #########################################
h3k4me3_late <- readRDS("./h3k4me3_late_2kbbin.rds")
label <- rownames(as.data.frame(h3k4me3_late@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k4me3_late@meta.data[["zga.signal.nor"]])

# 调整ZGA score
path <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/test/"
bed <- list.files(path)
k4_late_fc_titration <- list()
dif_portion <- c()
for (i in 1:length(bed)) {
  pathToBams1 <- c('./late2cell_h3k4me3/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
  regions <- paste0(path, bed[i])
  cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='late2cell_h3k4me3', min.cells = 0, min.regions = 0)
  h3k4me3_zga_late <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
  h3k4me3_zga_late <- colSums(h3k4me3_zga_late)
  
  cell_id <- read.csv("./late_k4me3_all_QC.txt", sep = "\t", header = F)
  cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
  h3k4me3_zga_late <- h3k4me3_zga_late[which(names(h3k4me3_zga_late) %in% cell_id$id)]
  h3k4me3_all_counts <- as.numeric(h3k4me3_late@meta.data[["nCount_peaks"]])
  label$zga_signal <- h3k4me3_zga_late/h3k4me3_all_counts
  
  embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
  k4_late_fc <- c()
  for (j in 1:length(embryos)) {
    dat_select <- label[grep(embryos[j], label$x5),]
    k4_late_fc_tmp <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
    k4_late_fc <- c(k4_late_fc, k4_late_fc_tmp)
  }
  k4_late_fc
  k4_late_fc_titration[[i]] <- k4_late_fc
  dif <- k4_late_fc[k4_late_fc > 1.3 | k4_late_fc < 0.8]
  dif_portion_tmp <- length(dif)/length(k4_late_fc)
  dif_portion <- c(dif_portion, dif_portion_tmp)
}
names(dif_portion) <- bed
dif_portion_k4_late <- dif_portion #majorZGA.genebody.5kb.bed, majorZGA.genebody.bed, MajorZGA_genes_tss5kb.bed, MajorZGA_genes_tss.bed均可达到52.8%


###### use MajorZGA_genes_tss.bed #####
label <- rownames(as.data.frame(h3k4me3_late@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k4me3_late@meta.data[["zga.signal.nor"]])

pathToBams1 <- c('./late2cell_h3k4me3/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/tss/MajorZGA_genes_tss.bed"
cisTopicObject_h3k4me3 <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='late2cell_h3k4me3', min.cells = 0, min.regions = 0)
h3k4me3_zga_late <- as.matrix(cisTopicObject_h3k4me3@count.matrix)
h3k4me3_zga_late <- colSums(h3k4me3_zga_late)

cell_id <- read.csv("./late_k4me3_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k4me3_zga_late <- h3k4me3_zga_late[which(names(h3k4me3_zga_late) %in% cell_id$id)]
h3k4me3_all_counts <- as.numeric(h3k4me3_late@meta.data[["nCount_peaks"]])
label$zga_signal <- h3k4me3_zga_late/h3k4me3_all_counts

embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
dim <- length(embryos) * 4
k4_late_zga <- matrix(1:dim, ncol = 4)
k4_late_zga <- as.data.frame(k4_late_zga)
for (j in 1:length(embryos)) {
  dat_select <- label[grep(embryos[j], label$x5),]
  k4_late_zga[j,1] <- embryos[j]
  k4_late_zga[j,2] <- dat_select$zga_signal[1]
  k4_late_zga[j,3] <- dat_select$zga_signal[2]
  k4_late_zga[j,4] <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
}
colnames(k4_late_zga) <- c("embryo", "cell1", "cell2", "foldChange")
k4_late_zga$low_zga <- pmin(k4_late_zga$cell1, k4_late_zga$cell2)
k4_late_zga$high_zga <- pmax(k4_late_zga$cell1, k4_late_zga$cell2)

library(ggplot2)
# 假设你的数据框是df，有ID列标识每行
df <- data.frame(ID = k4_late_zga$embryo, A = k4_late_zga$low_zga, B = k4_late_zga$high_zga)
# 将数据从宽格式转换为长格式
df_long <- reshape2::melt(df, id.vars = "ID", measure.vars = c("A", "B"))
pdf("./h3k4me3_late_c1_c2_zga_score.pdf")
ggplot(df_long, aes(x = variable, y = value, group = ID)) +
  geom_line(color = "lightgray") +  # 连接线
  geom_point() +  # 可选，增加点
  theme_bw() +  # 选择一个简洁的主题
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1.2) +  # 平均值连线
  stat_compare_means(comparisons = list(c("A", "B"))) +  # 统计显著性
  labs(x = "组别", y = "数值大小", title = "AB两列数值的连接线图")
dev.off()

  
  
  ################################### h3k27ac early #####################################
h3k27ac_early <- readRDS("./h3k27ac_early_5kbbin.rds")
label <- rownames(as.data.frame(h3k27ac_early@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k27ac_early@meta.data[["zga.signal.nor"]])

# 调整ZGA score
path <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/tss/"
bed <- list.files(path)
k27_early_fc_titration <- list()
dif_portion <- c()
for (i in 1:length(bed)) {
  pathToBams1 <- c('./early_h3k27ac/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
  regions <- paste0(path, bed[i])
  cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
  h3k27ac_zga_early <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
  h3k27ac_zga_early <- colSums(h3k27ac_zga_early)
  
  cell_id <- read.csv("./early_k27ac_all_QC.txt", sep = "\t", header = F)
  cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
  h3k27ac_zga_early <- h3k27ac_zga_early[which(names(h3k27ac_zga_early) %in% cell_id$id)]
  h3k27ac_all_counts <- as.numeric(h3k27ac_early@meta.data[["nCount_peaks"]])
  label$zga_signal <- h3k27ac_zga_early/h3k27ac_all_counts
  
  embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
  k27_early_fc <- c()
  for (j in 1:length(embryos)) {
    dat_select <- label[grep(embryos[j], label$x5),]
    k27_early_fc_tmp <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
    k27_early_fc <- c(k27_early_fc, k27_early_fc_tmp)
  }
  k27_early_fc
  k27_early_fc_titration[[i]] <- k27_early_fc
  dif <- k27_early_fc[k27_early_fc > 1.3 | k27_early_fc < 0.8]
  dif_portion_tmp <- length(dif)/length(k27_early_fc)
  dif_portion <- c(dif_portion, dif_portion_tmp)
}
names(dif_portion) <- bed
dif_portion_k27_early <- dif_portion
k27_early_fc_titration <- do.call(rbind, k27_early_fc_titration)

###### use MinorZGA_genes_tss.bed #####
label <- rownames(as.data.frame(h3k27ac_early@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k27ac_early@meta.data[["zga.signal.nor"]])

pathToBams1 <- c('./early_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/tss/MinorZGA_genes_tss.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='early2cell_h3k27ac', min.cells = 0, min.regions = 0)
h3k27ac_zga_early <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
h3k27ac_zga_early <- colSums(h3k27ac_zga_early)

cell_id <- read.csv("./early_k27ac_all_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k27ac_zga_early <- h3k27ac_zga_early[which(names(h3k27ac_zga_early) %in% cell_id$id)]
h3k27ac_all_counts <- as.numeric(h3k27ac_early@meta.data[["nCount_peaks"]])
label$zga_signal <- h3k27ac_zga_early/h3k27ac_all_counts

embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
dim <- length(embryos) * 4
k27_early_zga <- matrix(1:dim, ncol = 4)
k27_early_zga <- as.data.frame(k27_early_zga)
for (j in 1:length(embryos)) {
  dat_select <- label[grep(embryos[j], label$x5),]
  k27_early_zga[j,1] <- embryos[j]
  k27_early_zga[j,2] <- dat_select$zga_signal[1]
  k27_early_zga[j,3] <- dat_select$zga_signal[2]
  k27_early_zga[j,4] <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
}
colnames(k27_early_zga) <- c("embryo", "cell1", "cell2", "foldChange")
k27_early_zga$low_zga <- pmin(k27_early_zga$cell1, k27_early_zga$cell2)
k27_early_zga$high_zga <- pmax(k27_early_zga$cell1, k27_early_zga$cell2)

library(ggplot2)
# 假设你的数据框是df，有ID列标识每行
df <- data.frame(ID = k27_early_zga$embryo, A = k27_early_zga$low_zga, B = k27_early_zga$high_zga)
# 将数据从宽格式转换为长格式
df_long <- reshape2::melt(df, id.vars = "ID", measure.vars = c("A", "B"))
pdf("./h3k27me3_early_c1_c2_zga_score.pdf")
ggplot(df_long, aes(x = variable, y = value, group = ID)) +
  geom_line(color = "lightgray") +  # 连接线
  geom_point() +  # 可选，增加点
  theme_bw() +  # 选择一个简洁的主题
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1.2) +  # 平均值连线
  stat_compare_means(comparisons = list(c("A", "B"))) +  # 统计显著性
  labs(x = "组别", y = "数值大小", title = "AB两列数值的连接线图")
dev.off()

  
  ########################################### h3k27ac late #################################################
h3k27ac_late <- readRDS("./h3k27ac_late_5kbbin.rds")
label <- rownames(as.data.frame(h3k27ac_late@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k27ac_late@meta.data[["zga.signal.nor"]])

# 调整ZGA score
path <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/test/"
bed <- list.files(path)
k27_late_fc_titration <- list()
dif_portion <- c()
for (i in 1:length(bed)) {
  pathToBams1 <- c('./late2cell_h3k27ac/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
  regions <- paste0(path, bed[i])
  cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='late2cell_h3k27ac', min.cells = 0, min.regions = 0)
  h3k27ac_zga_late <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
  h3k27ac_zga_late <- colSums(h3k27ac_zga_late)
  
  cell_id <- read.csv("./late_k27acall_QC.txt", sep = "\t", header = F)
  cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
  h3k27ac_zga_late <- h3k27ac_zga_late[which(names(h3k27ac_zga_late) %in% cell_id$id)]
  h3k27ac_all_counts <- as.numeric(h3k27ac_late@meta.data[["nCount_peaks"]])
  label$zga_signal <- h3k27ac_zga_late/h3k27ac_all_counts
  
  embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
  k27_late_fc <- c()
  for (j in 1:length(embryos)) {
    dat_select <- label[grep(embryos[j], label$x5),]
    k27_late_fc_tmp <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
    k27_late_fc <- c(k27_late_fc, k27_late_fc_tmp)
  }
  k27_late_fc
  k27_late_fc_titration[[i]] <- k27_late_fc
  dif <- k27_late_fc[k27_late_fc > 1.3 | k27_late_fc < 0.8]
  dif_portion_tmp <- length(dif)/length(k27_late_fc)
  dif_portion <- c(dif_portion, dif_portion_tmp)
}
names(dif_portion) <- bed
dif_portion_k27_late <- dif_portion


###### use majorZGA.genebody.bed #####
label <- rownames(as.data.frame(h3k27ac_late@active.ident))
label <- as.data.frame(label)
label <- separate(label, label, c("x1","x2","x3","x4","x5","x6","x7", "x8"), sep = "_", remove = F)
label <- label[,c("label", "x5", "x6")]
label$zga_score <- as.numeric(h3k27ac_late@meta.data[["zga.signal.nor"]])

pathToBams1 <- c('./late2cell_h3k27ac/')
bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
bamFiles1 <- bamFiles1[grep("rmdup.bam", bamFiles1)]
regions <- "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/genebody/majorZGA.genebody.bed"
cisTopicObject_h3k27ac <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='late2cell_h3k27ac', min.cells = 0, min.regions = 0)
h3k27ac_zga_late <- as.matrix(cisTopicObject_h3k27ac@count.matrix)
h3k27ac_zga_late <- colSums(h3k27ac_zga_late)

cell_id <- read.csv("./late_k27acall_QC.txt", sep = "\t", header = F)
cell_id$id <- paste0(cell_id$V1, ".mm10_rmdup.bam")
h3k27ac_zga_late <- h3k27ac_zga_late[which(names(h3k27ac_zga_late) %in% cell_id$id)]
h3k27ac_all_counts <- as.numeric(h3k27ac_late@meta.data[["nCount_peaks"]])
label$zga_signal <- h3k27ac_zga_late/h3k27ac_all_counts

embryos <- unique(label$x5) # 统计同一个胚胎的两个细胞的zga score fold change.
dim <- length(embryos) * 4
k27_late_zga <- matrix(1:dim, ncol = 4)
k27_late_zga <- as.data.frame(k27_late_zga)
for (j in 1:length(embryos)) {
  dat_select <- label[grep(embryos[j], label$x5),]
  k27_late_zga[j,1] <- embryos[j]
  k27_late_zga[j,2] <- dat_select$zga_signal[1]
  k27_late_zga[j,3] <- dat_select$zga_signal[2]
  k27_late_zga[j,4] <- dat_select$zga_signal[1]/dat_select$zga_signal[2]
}
colnames(k27_late_zga) <- c("embryo", "cell1", "cell2", "foldChange")
k27_late_zga$low_zga <- pmin(k27_late_zga$cell1, k27_late_zga$cell2)
k27_late_zga$high_zga <- pmax(k27_late_zga$cell1, k27_late_zga$cell2)

library(ggplot2)
# 假设你的数据框是df，有ID列标识每行
df <- data.frame(ID = k27_late_zga$embryo, A = k27_late_zga$low_zga, B = k27_late_zga$high_zga)
# 将数据从宽格式转换为长格式
df_long <- reshape2::melt(df, id.vars = "ID", measure.vars = c("A", "B"))
pdf("./h3k27me3_late_c1_c2_zga_score.pdf")
ggplot(df_long, aes(x = variable, y = value, group = ID)) +
  geom_line(color = "lightgray") +  # 连接线
  geom_point() +  # 可选，增加点
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1.2) +  # 平均值连线
  stat_compare_means(comparisons = list(c("A", "B"))) +  # 统计显著性
  theme_bw() +  # 选择一个简洁的主题
  labs(x = "组别", y = "数值大小", title = "AB两列数值的连接线图")
dev.off()

write.table(k4_early_zga, "h3k4me3_early_zga_score.txt", sep = "\t", quote = F, row.names = F)
write.table(k4_late_zga, "h3k4me3_late_zga_score.txt", sep = "\t", quote = F, row.names = F)
write.table(k27_early_zga, "h3k27ac_early_zga_score.txt", sep = "\t", quote = F, row.names = F)
write.table(k27_late_zga, "h3k27ac_late_zga_score.txt", sep = "\t", quote = F, row.names = F)


# library(openxlsx)
# zga <- read.xlsx("/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/MajorZGA_genes.xlsx")
# zga <- zga$Gene
# zga <- na.omit(zga)
# ref <- read.csv("/media/helab/data1/min/00_reference/mm10_Refseq.transcript.uniq2.2.gtf", sep = "\t", header = F)
# ref$up <- ref$V4 -20000
# ref$down <- ref$V5 +20000
# zga_tss <- ref[which(ref$V9 %in% zga),]
# write.table(zga_tss[,c("V1", "up", "down")], "/media/helab/data1/min/00_reference/mouse_embryo/genes/ZGA/MajorZGA_genes_tss_20kb.bed",sep = "\t", quote = F, col.names = F, row.names = F)
