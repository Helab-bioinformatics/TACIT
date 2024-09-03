# setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/10_state_umap")
# library(tidyr)
# 
# mm10 <- read.csv("../08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)
# inte_1 <- read.csv("./chrhmm_200bin_loci_genebody.txt", sep = "\t", header = F)
# inte_2 <- read.csv("./chrhmm_200bin_loci_genebody_gene.txt", sep = "\t", header = F)
# inte_1$gene <- inte_2$V4
# inte_1 <- unite(inte_1, "loci", c("V1", "V2"), sep = ":", remove = F)
# inte_1 <- unite(inte_1, "loci2", c("loci", "V3"), sep = "-", remove = F)
# inte_1 <- inte_1[,-2]
# inte_1$loci <- mm10$V4[match(inte_1$loci2, mm10$V5)]
# write.table(inte_1, "mm10_genebody_chrhmm_loci.txt", quote = F, col.names = F, row.names = F, sep = "\t")  


################################# 4cell-1 ##########################
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_1/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E12", "E9")]
Pro <- mat[,c("E8")]
Te <- mat[,c("E10","E11")]
hetero_k9 <- mat[,c("E6","E4")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E5","E2")]

write.table(En, "./test/4cell_1_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_1_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_1_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_1_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_1_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_1_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E12", "E9")]
Pro <- mat[,c("E8")]
Te <- mat[,c("E10","E11")]
hetero_k9 <- mat[,c("E6","E4")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E5","E2")]

write.table(En, "./test/4cell_1_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_1_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_1_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_1_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_1_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_1_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


############################# 4cell-2 ############################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_2/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E1", "E3", "E4")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E11","E2")]
hetero_k9 <- mat[,c("E8")]
hetero_k27 <- mat[,c("E10")]
hetero <- mat[,c("E9","7")]

write.table(En, "./test/4cell_2_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_2_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_2_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_2_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_2_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_2_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E1", "E2","E3", "E4")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E11","E2")]
hetero_k9 <- mat[,c("E8")]
hetero_k27 <- mat[,c("E10")]
hetero <- mat[,c("E9","7")]

write.table(En, "./test/4cell_2_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_2_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_2_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_2_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_2_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_2_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


############################# 4cell-3 ############################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_3/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E12", "E11")]
Pro <- mat[,c("E1")]
Te <- mat[,c("E10","E9")]
hetero_k9 <- mat[,c("E6")]
hetero_k27 <- mat[,c("E4")]
hetero <- mat[,c("E3","E2")]

write.table(En, "./test/4cell_3_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_3_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_3_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_3_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_3_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_3_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E12", "E11")]
Pro <- mat[,c("E1")]
Te <- mat[,c("E10","E9")]
hetero_k9 <- mat[,c("E6")]
hetero_k27 <- mat[,c("E4")]
hetero <- mat[,c("E3","E2")]

write.table(En, "./test/4cell_3_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_3_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_3_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_3_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_3_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_3_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


############################### 4cell-4 ################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/03_chromHMM_4cell/4cell_4/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E1", "E4", "E2")]
Pro <- mat[,c("E3", "E5")]
Te <- mat[,c("E10")]
hetero_k9 <- mat[,c("E6")]
hetero_k27 <- mat[,c("E12")]
hetero <- mat[,c("E12","E6", "E7")]

write.table(En, "./test/4cell_4_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_4_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_4_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_4_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_4_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_4_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E1", "E4", "E2")]
Pro <- mat[,c("E3", "E5")]
Te <- mat[,c("E10")]
hetero_k9 <- mat[,c("E6")]
hetero_k27 <- mat[,c("E12")]
hetero <- mat[,c("E12","E6", "E7")]

write.table(En, "./test/4cell_4_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/4cell_4_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/4cell_4_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/4cell_4_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/4cell_4_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/4cell_4_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


###################### 2cell-1 ####################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/02_chromHMM_2cell/2cell_1/learnmodel_inte/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E12", "E11")]
Pro <- mat[,c("E9", "E1")]
Te <- mat[,c("E2", "E8")]
hetero_k9 <- mat[,c("E7")]
hetero_k27 <- mat[,c("E6", "E5")]
hetero <- mat[,c("E7")]

write.table(En, "./test/2cell_1_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/2cell_1_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/2cell_1_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/2cell_1_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/2cell_1_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/2cell_1_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E12", "E11")]
Pro <- mat[,c("E9", "E1")]
Te <- mat[,c("E2", "E8")]
hetero_k9 <- mat[,c("E7")]
hetero_k27 <- mat[,c("E6", "E5")]
hetero <- mat[,c("E7")]

write.table(En, "./test/2cell_1_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/2cell_1_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/2cell_1_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/2cell_1_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/2cell_1_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/2cell_1_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)

######################## 2cell_2 #################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/02_chromHMM_2cell/2cell_2/learnmodel_inte/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E12", "E9")]
Pro <- mat[,c("E1", "E11")]
Te <- mat[,c("E10", "E8")]
hetero_k9 <- mat[,c("E7")]
hetero_k27 <- mat[,c("E6", "E5")]
hetero <- mat[,c("E7")]

write.table(En, "./test/2cell_2_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/2cell_2_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/2cell_2_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/2cell_2_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/2cell_2_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/2cell_2_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E12", "E9")]
Pro <- mat[,c("E1", "E11")]
Te <- mat[,c("E10", "E8")]
hetero_k9 <- mat[,c("E7")]
hetero_k27 <- mat[,c("E6", "E5")]
hetero <- mat[,c("E7")]

write.table(En, "./test/2cell_2_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/2cell_2_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/2cell_2_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/2cell_2_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/2cell_2_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/2cell_2_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


########################### 8cell-1 (8cell鐨勫叏閮ㄧ敤8cell_all training鍚庣殑妯″瀷锛屾墍浠ラ兘鏄竴鏍风殑瀵瑰簲鍏崇郴)##########################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_1/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c( "E5", "E6")]
Pro <- mat[,c("E1", "E2")]
Te <- mat[,c("E7")]
hetero_k9 <- mat[,c("E12")]
hetero_k27 <- mat[,c("E9", "E10")]
hetero <- mat[,c("E11")]

write.table(En, "./test/8cell_1_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_1_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_1_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_1_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_1_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_1_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c( "E5", "E6")]
Pro <- mat[,c("E1", "E2")]
Te <- mat[,c("E7")]
hetero_k9 <- mat[,c("E12")]
hetero_k27 <- mat[,c("E9", "E10")]
hetero <- mat[,c("E11")]

write.table(En, "./test/8cell_1_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_1_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_1_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_1_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_1_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_1_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)

################################ 8cell_2 #################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_2/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E7", "E8", "E9", "E10")]
Pro <- mat[,c("E11", "E12")]
Te <- mat[,c("E6")]
hetero_k9 <- mat[,c("E1", "E2")]
hetero_k27 <- mat[,c("E4", "E2")]
hetero <- mat[,c("E1")]

write.table(En, "./test/8cell_2_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_2_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_2_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_2_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_2_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_2_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E7", "E8", "E9", "E10")]
Pro <- mat[,c("E11", "E12")]
Te <- mat[,c("E6")]
hetero_k9 <- mat[,c("E1", "E2")]
hetero_k27 <- mat[,c("E4", "E2")]
hetero <- mat[,c("E1")]

write.table(En, "./test/8cell_2_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_2_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_2_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_2_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_2_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_2_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


########################################### 8cell_3 ###################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_3/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E4", "E6", "E7", "E8")]
Pro <- mat[,c("E1")]
Te <- mat[,c("E3", "E5")]
hetero_k9 <- mat[,c("E12", "E9", "E11")]
hetero_k27 <- mat[,c("E9")]
hetero <- mat[,c("E12")]

write.table(En, "./test/8cell_3_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_3_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_3_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_3_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_3_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_3_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E4", "E6", "E7", "E8")]
Pro <- mat[,c("E1")]
Te <- mat[,c("E3", "E5")]
hetero_k9 <- mat[,c("E12","E9", "E11")]
hetero_k27 <- mat[,c("E9")]
hetero <- mat[,c("E12")]

write.table(En, "./test/8cell_3_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_3_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_3_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_3_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_3_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_3_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


############################### 8cell_4 ##########################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_4/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E8", "E10", "E11")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E7", "E9")]
hetero_k9 <- mat[,c("E3", "E4", "E5")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E1")]

write.table(En, "./test/8cell_4_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_4_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_4_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_4_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_4_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_4_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E8", "E10", "E11")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E7", "E9")]
hetero_k9 <- mat[,c("E3", "E4", "E5")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E1")]

write.table(En, "./test/8cell_4_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_4_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_4_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_4_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_4_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_4_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


########################################### 8cell_5 #################################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_5/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E7", "E10", "E11")]
Pro <- mat[,c("E9", "E12")]
Te <- mat[,c("E6","E8")]
hetero_k9 <- mat[,c("E4")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E3")]

write.table(En, "./test/8cell_5_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_5_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_5_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_5_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_5_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_5_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E7", "E10", "E11")]
Pro <- mat[,c("E9", "E12")]
Te <- mat[,c("E6","E8")]
hetero_k9 <- mat[,c("E4")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E3")]

write.table(En, "./test/8cell_5_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_5_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_5_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_5_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_5_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_5_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


########################################### 8cell_6 ##############################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_6/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E9", "E10", "E11")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E7")]
hetero_k9 <- mat[,c("E4")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E2", "E5")]

write.table(En, "./test/8cell_6_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_6_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_6_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_6_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_6_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_6_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E9", "E10", "E11")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E7")]
hetero_k9 <- mat[,c("E4")]
hetero_k27 <- mat[,c("E1")]
hetero <- mat[,c("E2", "E5")]

write.table(En, "./test/8cell_6_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_6_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_6_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_6_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_6_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_6_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


################################### 8cell_7 #############################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_7/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E2", "E5")]
Pro <- mat[,c("E1")]
Te <- mat[,c("E3", "E4", "E6", "E7", "E9")]
hetero_k9 <- mat[,c("E11")]
hetero_k27 <- mat[,c("E12")]
hetero <- mat[,c("E12")]

write.table(En, "./test/8cell_7_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_7_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_7_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_7_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_7_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_7_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E2", "E5")]
Pro <- mat[,c("E1")]
Te <- mat[,c("E3", "E4", "E6", "E7", "E9")]
hetero_k9 <- mat[,c("E11")]
hetero_k27 <- mat[,c("E12")]
hetero <- mat[,c("E12")]

write.table(En, "./test/8cell_7_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_7_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_7_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_7_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_7_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_7_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


########################################## 8cell_8 #######################################################
rm(list = ls())
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/04_chromHMM_8cell/8cell_8/learnmodel/POSTERIOR")
mat <- read.csv("./test/all_posterior.txt", sep = "\t")
En <- mat[,c("E1", "E2", "E3")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E4", "E5")]
hetero_k9 <- mat[,c("E10", "E11")]
hetero_k27 <- mat[,c("E8")]
hetero <- mat[,c("E9")]

write.table(En, "./test/8cell_8_enhancer_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_8_promoter_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_8_transcription_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_8_heterochromatin_h3k9me3_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_8_heterochromatin_h3k27me3_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_8_heterochromatin_both_posterior.txt", sep="\t", quote = F)

mat <- read.csv("./test/all_posterior.txt", sep = "\t")
tss <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/mm10_tss2kb_chrhmm_loci.txt", sep = "\t", header = F)
mat <- mat[tss$V6,]
En <- mat[,c("E1", "E2", "E3")]
Pro <- mat[,c("E12")]
Te <- mat[,c("E4", "E5")]
hetero_k9 <- mat[,c("E10", "E11")]
hetero_k27 <- mat[,c("E8")]
hetero <- mat[,c("E9")]

write.table(En, "./test/8cell_8_enhancer_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Pro, "./test/8cell_8_promoter_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(Te, "./test/8cell_8_transcription_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k9, "./test/8cell_8_heterochromatin_h3k9me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero_k27, "./test/8cell_8_heterochromatin_h3k27me3_tss2kb_posterior.txt", sep="\t", quote = F)
write.table(hetero, "./test/8cell_8_heterochromatin_tss2kb_posterior.txt", sep="\t", quote = F)


