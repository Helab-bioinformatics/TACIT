setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/02_tss_umap/")

############################## heterochromatin_h3k9me3 strong ################################
rm(list = ls())
zygote <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/01_chromHMM_zygote/learnmodel_nochrM_post/POSTERIOR/test/zygote_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C2_1 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/02_chromHMM_2cell/2cell_1/learnmodel_inte/POSTERIOR/test/2cell_1_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C2_2 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/02_chromHMM_2cell/2cell_2/learnmodel_inte/POSTERIOR/test/2cell_2_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C4_1 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_1/learnmodel/POSTERIOR/test/4cell_1_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C4_2 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_2/learnmodel/POSTERIOR/test/4cell_2_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C4_3 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_3/learnmodel/POSTERIOR/test/4cell_3_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C4_4 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_4/learnmodel/POSTERIOR/test/4cell_4_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_1 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_1/learnmodel/POSTERIOR/test/8cell_1_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_2 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_2/learnmodel/POSTERIOR/test/8cell_2_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_3 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_3/learnmodel/POSTERIOR/test/8cell_3_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_4 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_4/learnmodel/POSTERIOR/test/8cell_4_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_5 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_5/learnmodel/POSTERIOR/test/8cell_5_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_6 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_6/learnmodel/POSTERIOR/test/8cell_6_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_7 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_7/learnmodel/POSTERIOR/test/8cell_7_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")
C8_8 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_8/learnmodel/POSTERIOR/test/8cell_8_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t")

zygote <- apply(zygote,1,mean)
C2_1 <- apply(C2_1,1,mean)
C2_2 <- apply(C2_2,1,mean)
C4_1 <- apply(C4_1,1,mean)
C4_2 <- apply(C4_2,1,mean)
C4_3 <- apply(C4_3,1,mean)
C4_4 <- apply(C4_4,1,mean)
C8_1 <- apply(C8_1,1,mean)
C8_2 <- apply(C8_2,1,mean)
C8_3 <- apply(C8_3,1,mean)
C8_4 <- apply(C8_4,1,mean)
C8_5 <- apply(C8_5,1,mean)
C8_6 <- apply(C8_6,1,mean)
C8_7 <- apply(C8_7,1,mean)
C8_8 <- apply(C8_8,1,mean)
mat <- data.frame(zygote=zygote, C2_1=C2_1, C2_2=C2_2,
                  C4_1=C4_1, C4_2=C4_2, C4_3=C4_3, C4_4=C4_4,
                  C8_1=C8_1, C8_2=C8_2,C8_3=C8_3, C8_4=C8_4,C8_5=C8_5, C8_6=C8_6,C8_7=C8_7, C8_8=C8_8)
rm(zygote, C2_1, C2_2, C4_1, C4_2, C4_3, C4_4, C8_1, C8_2, C8_3, C8_4, C8_5, C8_6, C8_7, C8_8)
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
nrow(mat)
nrow(tss)
mat$loci <- tss$V6
mat <- mat[!duplicated(mat$loci),]
rownames(mat) <- mat$loci
write.table(mat[,1:15], "all_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t", quote = F)

############################ 合并一样的行，且去掉单独的行 #########################
library(matrixStats)
mat <- read.csv("all_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep = "\t" )
mat <- mat[,-1]
mat <- mat[which(rowSums(mat) > 0),] #去掉所有细胞中都是概率为0的位置
mat_rou <- round(mat, 2)
rownames(mat_rou) <- rownames(mat)
mat_rou <- mat_rou[which(rowSums(mat_rou) > 0),]
mat_rou <- as.matrix(mat_rou)
mat_rou <- mat_rou[which(rowSds(mat_rou) >0),]
mat_rou <- as.data.frame(mat_rou)

# 创建一个逻辑向量，表示当前行与前一行的值是否完全相等
is_equal_to_previous <- c(FALSE, apply(mat_rou[-1, ] == mat_rou[-nrow(mat_rou), ], 1, all))
# 初始化deleted_rows向量
deleted_rows <- numeric(length = nrow(mat_rou))
# 遍历数据框，计算每个相邻范围内相同的行数（删除的行数）
start_range <- 1
for (i in 2:nrow(mat_rou)) {
  if (!is_equal_to_previous[i]) {
    deleted_rows[start_range] <- sum(is_equal_to_previous[start_range:(i - 1)])
    start_range <- i
  }
}
# 只保留与前一行不完全相等的行
filtered_mat_rou <- mat_rou[!is_equal_to_previous, ]
# 添加记录删除行数的第16列
filtered_mat_rou$deleted_rows <- deleted_rows[!is_equal_to_previous]
# 保留原始的行名
rownames(filtered_mat_rou) <- rownames(mat_rou)[!is_equal_to_previous]

filtered_mat_rou <- filtered_mat_rou[filtered_mat_rou$deleted_rows != 0,]
filtered_mat_rou <- filtered_mat_rou[filtered_mat_rou$deleted_rows != 1,]

write.table(filtered_mat_rou, "all_heterochromatin_h3k9me3_tss2kb_posterior_selected.txt", sep = "\t", quote = F)





