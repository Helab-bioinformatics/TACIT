setwd("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/03_random_forest")

zygote <- read.csv("./zygote_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C2_1 <- read.csv("./2cell_1_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C2_2 <- read.csv("./2cell_2_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C4_1 <- read.csv("./4cell_1_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C4_2 <- read.csv("./4cell_2_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C4_3 <- read.csv("./4cell_3_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C4_4 <- read.csv("./4cell_4_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_1 <- read.csv("./8cell_1_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_2 <- read.csv("./8cell_2_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_3 <- read.csv("./8cell_3_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_4 <- read.csv("./8cell_4_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_5 <- read.csv("./8cell_5_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_6 <- read.csv("./8cell_6_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_7 <- read.csv("./8cell_7_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)
C8_8 <- read.csv("./8cell_8_mm10_chrhmm_dense_anno_200.txt", sep = "\t", header = F)

library(tidyr)
cells <- list(zygote, C2_1, C2_2, C4_1, C4_2,C4_3, C4_4, C8_1, C8_2,C8_3, C8_4,C8_5, C8_6,C8_7, C8_8)
for (i in 1:length(cells)) {
  a <- cells[[i]]
  a <- unite(a, loci, c("V1", "V2", "V3"), sep = "-")
  cells[[i]] <- a
}

inte <- Reduce(intersect, list(cells[[1]]$loci, cells[[2]]$loci, cells[[3]]$loci, cells[[4]]$loci, cells[[5]]$loci, cells[[6]]$loci, cells[[7]]$loci, cells[[8]]$loci, 
                               cells[[9]]$loci, cells[[10]]$loci, cells[[11]]$loci, cells[[12]]$loci, cells[[13]]$loci, cells[[14]]$loci, cells[[15]]$loci))
for (i in 1:length(cells)) {
  a <- cells[[i]]
  a <- a[which(a$loci %in% inte),]
  cells[[i]] <- a
}

anno <- data.frame(zygote=cells[[1]]$V4, C2_1=cells[[2]]$V4, C2_2=cells[[3]]$V4,
                   C4_1=cells[[4]]$V4, C4_2=cells[[5]]$V4, C4_3=cells[[6]]$V4, C4_4=cells[[7]]$V4,
                   C8_1=cells[[8]]$V4, C8_2=cells[[9]]$V4, C8_3=cells[[10]]$V4, C8_4=cells[[11]]$V4, C8_5=cells[[12]]$V4, C8_6=cells[[13]]$V4, C8_7=cells[[14]]$V4, C8_8=cells[[15]]$V4)
rownames(anno) <- inte
anno <- anno[which(rowSums(anno == "Un") <13),]
write.table(anno, "./zygoteTo8cell_all_anno_2.txt", sep = "\t", quote = F) 

anno <- read.csv("./zygoteTo8cell_all_anno.txt", sep = "\t")
############################################
is_equal_to_previous <- c(FALSE, apply(anno[-1, ] == anno[-nrow(anno), ], 1, all))
deleted_rows <- numeric(length = nrow(anno))
start_range <- 1
for (i in 2:nrow(anno)) {
  if (!is_equal_to_previous[i]) {
    deleted_rows[start_range] <- sum(is_equal_to_previous[start_range:(i - 1)])
    start_range <- i
  }
}
filtered_anno <- anno[!is_equal_to_previous, ]
filtered_anno$deleted_rows <- deleted_rows[!is_equal_to_previous]
rownames(filtered_anno) <- rownames(anno)[!is_equal_to_previous]

filtered_anno <- filtered_anno[which(filtered_anno$deleted_rows >2),]
filtered_anno <- filtered_anno[which(rowSums(filtered_anno[,1:6] == "Un") <5),]
write.table(filtered_anno, "zygoteTo8cell_filtered_anno_2.txt", sep = "\t", quote = F)


mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/08_repressive/chrhmm_200bin_loci.txt", sep = "\t", header = F)
mm10 <- unite(mm10, loci, c("V1", "V2", "V3"), sep = "-")
filtered_anno$loci <- mm10$V4[match(rownames(filtered_anno), mm10$loci)]
filtered_anno <- filtered_anno[!duplicated(filtered_anno$loci),]
write.table(filtered_anno, "zygoteTo8cell_filtered_anno_2.txt", sep = "\t", quote = F)

# filtered_anno <- filtered_anno[which(filtered_anno$deleted_rows <300),]
# filtered_anno <- filtered_anno[which(rowSums(filtered_anno[,1:2] == "Un") <2),]

############################## RNA ##################################
RNA <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/RNA_pseudo.rds")
genes <- as.data.frame(RNA@assays[["RNA"]]@counts)
toti <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/totipotency.markers.csv", sep = " ", header = F)

genes_toti <- genes[which(rownames(genes) %in% toti$V1),]
genes_toti <- genes_toti[,grep("zygote|2cell|4cell|8cell", colnames(genes_toti))]
genes_toti <- t(as.matrix(genes_toti))
genes_toti <- scale(genes_toti)
genes_toti <- as.data.frame(t(genes_toti))

genes_toti <- na.omit(genes_toti)
genes_zygote <- genes_toti[,grep("zygote", colnames(genes_toti))]
genes_2cellC1 <- genes_toti[,c("2cell_pseudo_1", "2cell_pseudo_6", "2cell_pseudo_8", "2cell_pseudo_9", "2cell_pseudo_13")]
genes_2cellC2 <- genes_toti[,c("2cell_pseudo_14", "2cell_pseudo_22", "2cell_pseudo_24", "2cell_pseudo_25", "2cell_pseudo_26")]
genes_4cellC1 <- genes_toti[,c("4cell_pseudo_1", "4cell_pseudo_2", "4cell_pseudo_3", "4cell_pseudo_5", "4cell_pseudo_6")]
genes_4cellC2 <- genes_toti[,c("4cell_pseudo_9", "4cell_pseudo_11", "4cell_pseudo_13", "4cell_pseudo_14")]
genes_4cellC3 <- genes_toti[,c("4cell_pseudo_15", "4cell_pseudo_16", "4cell_pseudo_17", "4cell_pseudo_18")]
genes_4cellC4 <- genes_toti[,c("4cell_pseudo_19", "4cell_pseudo_20", "4cell_pseudo_22", "4cell_pseudo_23")]
genes_8cell <- genes_toti[,grep("8cell", colnames(genes_toti))]

genes_2cellC1 <- rowMeans(genes_2cellC1)
genes_2cellC2 <- rowMeans(genes_2cellC2)
genes_4cellC1 <- rowMeans(genes_4cellC1)
genes_4cellC2 <- rowMeans(genes_4cellC2)
genes_4cellC3 <- rowMeans(genes_4cellC3)
genes_4cellC4 <- rowMeans(genes_4cellC4)
genes_8cell <- rowMeans(genes_8cell)

exp <- c(median(genes_2cellC1), median(genes_2cellC2), median(genes_4cellC1), median(genes_4cellC2),
         median(genes_4cellC3), median(genes_4cellC4), median(genes_8cell))


############################# pearson correlation between genes and chromatin states ####################################
filtered_anno <- read.csv("./zygoteTo8cell_filtered_anno.txt", sep = "\t")

promoter <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/01_allbin_umap/all_promoter_posterior_selected.txt", sep = "\t")
ens <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/01_allbin_umap/all_enhancer_strong_posterior_selected.txt", sep = "\t")
tran <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/01_allbin_umap/all_transcription_posterior_selected.txt", sep = "\t")
k27 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/01_allbin_umap/all_heterochromatin_h3k27me3_posterior_selected.txt", sep = "\t")
k9 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/07_fig_3/01_allbin_umap/all_heterochromatin_h3k9me3_posterior_selected.txt", sep = "\t")


library(tidyr)
################### promoter ####################
promoter <- promoter[which(rownames(promoter) %in% filtered_anno$loci),]
promoter$C8_all <- apply(promoter[,7:14], 1, median)
promoter <- promoter[,c(1:6, 16)]

promoter <- as.matrix(promoter)
library(stats)
cor <- c()
for (i in 1:nrow(promoter)) {
  cor[i] <- cor(exp, promoter[i,], method = 'pearson')
}
names(cor) <- rownames(promoter)
cor_promoter <- na.omit(cor)
cor_promoter_high <- cor_promoter[cor_promoter > 0.3]
cor_promoter_low <- cor_promoter[cor_promoter < -0.3]
cor_promoter_no <- cor_promoter[cor_promoter >= -0.3 & cor_promoter <= 0.5]

#write.table(cor_promoter, "promoter_filtered_pearson_correlation.txt", sep = "\t", quote = F)

################### enhancer poised #######################
ens <- ens[which(rownames(ens) %in% filtered_anno$loci),]
ens$C8_all <- apply(ens[,7:14], 1, median)
ens <- ens[,c(1:6, 16)]

ens <- as.matrix(ens)
library(stats)
cor <- c()
for (i in 1:nrow(ens)) {
  cor[i] <- cor(exp, ens[i,], method = 'pearson')
}
names(cor) <- rownames(ens)
cor_ens <- na.omit(cor)
cor_ens_high <- cor_ens[cor_ens > 0.3]
cor_ens_low <- cor_ens[cor_ens < -0.3]
cor_ens_no <- cor_ens[cor_ens >= -0.3 & cor_ens <= 0.3]

#write.table(cor_enp, "enp_filtered_pearson_correlation.txt", sep = "\t", quote = F)


################### transcription #######################
tran <- tran[which(rownames(tran) %in% filtered_anno$loci),]
tran$C8_all <- apply(tran[,7:14], 1, median)
tran <- tran[,c(1:6, 16)]

tran <- as.matrix(tran)
library(stats)
cor <- c()
for (i in 1:nrow(tran)) {
  cor[i] <- cor(exp, tran[i,], method = 'pearson')
}
names(cor) <- rownames(tran)
cor_tran <- na.omit(cor)
cor_tran_high <- cor_tran[cor_tran > 0.3]
cor_tran_low <- cor_tran[cor_tran < -0.3]
cor_tran_no <- cor_tran[cor_tran >= -0.3 & cor_tran <= 0.3]
#write.table(cor_tran, "tran_filtered_pearson_correlation.txt", sep = "\t", quote = F)


################### heterochromatin #######################
k9 <- k9[which(rownames(k9) %in% filtered_anno$loci),]##无overlap
k9$C8_all <- apply(k9[,7:14], 1, median)
k9 <- k9[,c(1:6, 16)]

k9 <- as.matrix(k9)
library(stats)
cor <- c()
for (i in 1:nrow(k9)) {
  cor[i] <- cor(exp, k9[i,], method = 'pearson')
}
names(cor) <- rownames(k9)
cor_k9 <- na.omit(cor)
cor_k9_high <- cor_k9[cor_k9 > 0.3]
cor_k9_low <- cor_k9[cor_k9 < -0.3]
cor_k9_no <- cor_k9[cor_k9 >= -0.3 & cor_k9 <= 0.3]

k27 <- k27[which(rownames(k27) %in% filtered_anno$loci),]##无overlap
k27$C8_all <- apply(k27[,7:14], 1, median)
k27 <- k27[,c(1:6, 16)]

k27 <- as.matrix(k27)
library(stats)
cor <- c()
for (i in 1:nrow(k27)) {
  cor[i] <- cor(exp, k27[i,], method = 'pearson')
}
names(cor) <- rownames(k27)
cor_k27 <- na.omit(cor)
cor_k27_high <- cor_k27[cor_k27 > 0.3]
cor_k27_low <- cor_k27[cor_k27 < -0.3]
cor_k27_no <- cor_k27[cor_k27 >= -0.3 & cor_k27 <= 0.3]

#######################################################################
cor <- c(cor_promoter_high, cor_promoter_low, cor_ens_high, cor_ens_low, cor_tran_high, cor_tran_low,
         cor_k27_high, cor_k27_low, cor_k9_high, cor_k9_low)
loci <- names(cor)
loci <- unique(loci)

cor_no <- c(cor_promoter_no, cor_ens_no, cor_tran_no, cor_k27_no, cor_k9_no)

mat <- filtered_anno[which(filtered_anno$loci %in% loci),1:15]
rownames(mat) <- filtered_anno[which(filtered_anno$loci %in% loci),17]
write.table(mat, "zygoteTo8cell_filtered_high_corredlated.txt", sep="\t", quote=F)

##########################这里用了第一次挑出的6331 loci进行traning，但是用的是第二次注释的结果，即2cell_1的E9改成Ts ################
filtered_anno <- read.csv("./zygoteTo8cell_filtered_anno_2.txt", sep = "\t")
mat <- read.csv("./zygoteTo8cell_filtered_high_corredlated.txt", sep = "\t")
mat <- filtered_anno[which(filtered_anno$loci %in% rownames(mat)),]
rownames(mat) <- mat$loci
mat <- mat[,1:15]

#mat <- read.csv("./zygoteTo8cell_filtered_high_corredlated.txt", sep = "\t")
mat <- na.omit(mat)
mat[mat=="Ps"] <- 1
mat[mat=="Pw"] <- 1
mat[mat=="Es"] <- 2
mat[mat=="Ew"] <- 2
mat[mat=="Tw"] <- 3
mat[mat=="Ts"] <- 3
mat[mat=="K27"] <- 4
mat[mat=="K9"] <- 4
mat[mat=="He"] <- 4
mat[mat=="Un"] <- 5
mat[mat=="Mu"] <- 5

mat <- mat[,-1]
row <- rownames(mat)
mat <- apply(mat, 2, as.numeric)
rownames(mat) <- row
mat <- as.matrix(mat)
toti <- mat[,c("C2_1","C2_2","C8_1","C8_2","C8_3","C8_4","C8_5","C8_6","C8_7","C8_8")]

################################# randomForest ###############################
library(randomForest)
library(pROC)
library(ROCR)
library(matrixStats)
library(tidyr)
set.seed(12)

toti <- as.data.frame(t(toti))
toti$identity <- as.factor(c(rep("high", 2), rep("no", 8)))

train <- toti
err<-as.numeric()
set.seed(12)
#options(expressions = 5e5)
for(i in c(1:10)){
  mtry_test <- randomForest(identity ~., data=train,ntree=10000,mtry=i)
  err<- append( err, mean( mtry_test$err.rate ) )
}
mtry<-which.min(err)#mtry=3
set.seed(123)
mtry_fit <- randomForest(identity ~.,data = train, mtry=mtry, ntree=10000, importance=T, proximity=TRUE)
mtry_fit #mtry=3, ntree=10000, OOB error=0%


pdf("./train_all_filtered_states.pdf")
plot(mtry_fit)#best ntree=10000
plot(mtry_fit$err.rate)
title(main = "use all filtered bins", sub = "err=10%")

#
importance_otu <- importance(mtry_fit)
varImpPlot(mtry_fit)
importance_otu <- as.data.frame(importance_otu)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)
write.table(importance_otu, './all_importance_otu_2.txt', sep = '\t', quote = FALSE)
dev.off()
#write.table(importance_otu[1:30, ], 'importance_otu_top30.txt', sep = '\t', col.names = NA, quote = FALSE)
# 

set.seed(123)
otu_train.cv <- replicate(5, rfcv(train[-ncol(train)], train$identity, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

library(ggplot2)
ibrary(splines) 
pdf("test.pdf")
p <- ggplot(otu_train.cv, aes(otus, value)) +
   geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
   labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
p
dev.off()

saveRDS(mtry_fit, "./toti_random_2.rds")
write.table(toti, "./toti_2cell_8cell_training.txt", sep = "\t", quote = F)

#################################### heatmap for 4200 bins ####################################################
library(matrixStats)
library(tidyr)
toti <- read.csv("./toti_2cell_8cell_training.txt", sep = "\t")
importance_otu <- read.csv("./all_importance_otu.txt", sep = "\t")
toti <- toti[,rownames(importance_otu)[1:2400]]
toti <- t(as.matrix(toti))

set.seed(12)
toti <- toti[which(rowSds(toti)>0),]
km<- kmeans(toti,10) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(toti, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
m.kmeans<- m.kmeans[o,]# order the matrix according to the order of the last column
anno_row <- as.data.frame(m.kmeans[,11])
anno_row[,2] <- anno_row[,1]

library(pheatmap)
pdf("top2400_kmeans_heatmap.pdf")
pheatmap( m.kmeans[,1:10],
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.05, annotation_row = anno_row[,1:2],
          show_rownames=FALSE, show_colnames=FALSE,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

m.kmeans <- as.data.frame(m.kmeans)
m.kmeans_1 <- m.kmeans[which(m.kmeans$V11 == "1"),]
m.kmeans_2 <- m.kmeans[grep("2", m.kmeans$V11),]
m.kmeans_3 <- m.kmeans[grep("3", m.kmeans$V11),]
m.kmeans_4 <- m.kmeans[grep("4", m.kmeans$V11),]
m.kmeans_5 <- m.kmeans[grep("5", m.kmeans$V11),]
m.kmeans_6 <- m.kmeans[grep("6", m.kmeans$V11),]
m.kmeans_7 <- m.kmeans[grep("7", m.kmeans$V11),]
m.kmeans_8 <- m.kmeans[grep("8", m.kmeans$V11),]
m.kmeans_9 <- m.kmeans[grep("9", m.kmeans$V11),]
m.kmeans_10 <- m.kmeans[grep("10", m.kmeans$V11),]
m.kmeans_re <- rbind(m.kmeans_6, m.kmeans_9)
m.kmeans_re <- m.kmeans_re[sample(rownames(m.kmeans_re)),]
m.kmeans_new_order <- rbind(m.kmeans_8, m.kmeans_7, m.kmeans_2, m.kmeans_3, 
                            m.kmeans_10, m.kmeans_5, m.kmeans_1, m.kmeans_4, m.kmeans_re)
pdf("top1920_kmeans_heatmap_2_new_order.pdf")
pheatmap( m.kmeans_new_order[,1:10],
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.05, #annotation_row = anno_row[,1:2],
          show_rownames=FALSE, show_colnames=FALSE,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

write.table(rownames(m.kmeans_new_order), "top1920_kmeans_order.txt", sep = "\t", quote = F)

################## calculate similarity between 2-cell and 4-cell ##################
filtered_anno <- read.csv("./zygoteTo8cell_filtered_anno_2.txt", sep = "\t")
mat <- read.csv("./zygoteTo8cell_filtered_high_corredlated.txt", sep = "\t")
toti <- read.csv("./toti_2cell_8cell_training.txt", sep = "\t")
toti <- t(as.matrix(toti))
toti <- toti[rownames(m.kmeans_new_order),]
mat <- mat[rownames(m.kmeans_new_order),]
mat_new <- cbind(mat[,1], toti[,1:2], mat[,4:15])

mat_new <- na.omit(mat_new)
mat_new[mat_new=="Ps"] <- 1
mat_new[mat_new=="Pw"] <- 1
mat_new[mat_new=="Es"] <- 2
mat_new[mat_new=="Ew"] <- 2
mat_new[mat_new=="Tw"] <- 3
mat_new[mat_new=="Ts"] <- 3
mat_new[mat_new=="K27"] <- 4
mat_new[mat_new=="K9"] <- 4
mat_new[mat_new=="He"] <- 4
mat_new[mat_new=="Un"] <- 5
mat_new[mat_new=="Mu"] <- 5

row <- rownames(mat_new)
mat_new <- apply(mat_new, 2, as.numeric)
rownames(mat_new) <- row
mat_new <- as.matrix(mat_new)
pdf("top1920_kmeans_heatmap_2_new_order_allcells.pdf")
pheatmap( mat_new,
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.05, #annotation_row = anno_row[,1:2],
          show_rownames=FALSE, show_colnames=FALSE,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()


cramerV = function(x, y=NULL, 
                   ci=FALSE, conf=0.95, type="perc",
                   R=1000, histogram=FALSE, 
                   digits=4, bias.correct=FALSE, 
                   reportIncomplete=FALSE, 
                   verbose=FALSE, ...) {
  
  CV=NULL
  
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){y=as.vector(y)}
  if(is.vector(x) & is.vector(y)){
    N      = length(x)
    Chi.sq = suppressWarnings(chisq.test(x, y, correct=FALSE, ...)$statistic)
    Phi    = Chi.sq / N
    Row    = length(unique(x))
    C      = length(unique(y))
    CV     =  sqrt(Phi / min(Row-1, C-1))
  }
  
  if(is.matrix(x)){x=as.table(x)}
  if(is.table(x)){
    TABLE  = x
    N      = sum(TABLE)
    Chi.sq = suppressWarnings(chisq.test(TABLE, correct=FALSE, ...)$statistic)
    Phi    = Chi.sq / N
    Row    = nrow(x)
    C      = ncol(x)
    CV     =  sqrt(Phi / min(Row-1, C-1))
  }
  
  PhiOrg = Phi
  VOrg   = CV
  PhiNew = NA
  VNew   = NA
  if(bias.correct){
    Phi     = max(0, Phi-((Row-1)*(C-1)/(N-1)))
    CC      = C-((C-1)^2/(N-1))
    RR     = Row-((Row-1)^2/(N-1))
    CV     = sqrt(Phi / min(RR-1, CC-1))
    PhiNew = Phi
    VNew   = CV
  } 
  
  if(verbose){
    cat("\n")
    cat("Rows          =", signif(Row, digits=digits))
    cat("\n")
    cat("Columns       =", signif(C, digits=digits))
    cat("\n")
    cat("N             =", signif(N, digits=digits))
    cat("\n")
    cat("Chi-squared   =", signif(Chi.sq, digits=digits))
    cat("\n")
    cat("Phi           =", signif(PhiOrg, digits=digits))
    cat("\n")
    cat("Corrected Phi =", signif(PhiNew, digits=digits))
    cat("\n")
    cat("V             =", signif(VOrg, digits=digits))
    cat("\n")
    cat("Corrected V   =", signif(VNew, digits=digits))
    cat("\n")     
    cat("\n")          
  }
  
  if(bias.correct){
    PhiNew = max(0, Phi-((Row-1)*(C-1)/(N-1)))
    CC  = C-((C-1)^2/(N-1))
    RR  = Row-((Row-1)^2/(N-1))
    CV  = sqrt(Phi / min(RR-1, CC-1))
  }
  
  CV = signif(as.numeric(CV), digits=digits)
  
  if(is.nan(CV) & ci==TRUE){
    return(data.frame(Cramer.V=CV, lower.ci=NA, upper.ci=NA))} 
  
  if(ci==TRUE){
    if(is.matrix(x)){x=as.table(x)}
    if(is.table(x)){
      Counts = as.data.frame(x)
      Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
      rownames(Long) = seq(1:nrow(Long))
    }
    if(is.vector(x) & is.vector(y)){  
      Long = data.frame(x=x, y=y)
    }
    
    L1     = length(unique(droplevels(Long[,1])))
    L2     = length(unique(droplevels(Long[,2])))
    
    Function = function(input, index){
      Input = input[index,]
      
      NOTEQUAL=0
      if(length(unique(droplevels(Input[,1]))) != L1 |
         length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}
      
      if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
      
      if(NOTEQUAL==0){
        N      = length(Input[,1])
        Chi.sq = suppressWarnings(chisq.test(Input[,1], Input[,2], 
                                             correct=FALSE, ...)$statistic)
        Phi    =  Chi.sq / N
        Row    = length(unique(Input[,1]))
        C    = length(unique(Input[,2]))
        CV     =  sqrt(Phi / min(Row-1, C-1))
        FLAG   =  0
        
        if(bias.correct==TRUE){
          Phi = max(0, Phi-((Row-1)*(C-1)/(N-1)))
          CC  = C-((C-1)^2/(N-1))
          RR  = Row-((Row-1)^2/(N-1))
          CV  = sqrt(Phi / min(RR-1, CC-1))
        }
        
        return(c(CV,FLAG))
      }
    }
    
    Boot = boot(Long, Function, R=R)
    BCI  = boot.ci(Boot, conf=conf, type=type)
    if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
    if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
    if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
    if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}
    
    if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
    
    CI1=signif(CI1, digits=digits)
    CI2=signif(CI2, digits=digits)
    
    if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="V", main="")}
    
  }
  if(ci==FALSE){names(CV)="Cramer V"; return(CV)}
  if(ci==TRUE){return(data.frame(Cramer.V=CV, lower.ci=CI1, upper.ci=CI2))}  
}

# 
cramer_matrix <- matrix(1:225,ncol = 15)
for (i in 1:15) {
  for (j in i:15) {
    cramer_matrix[i,j] <- cramerV(mat_new_c8[,i], mat_new_c8[,j])
    cramer_matrix[j,i] <- cramer_matrix[i,j]
  }
}
rownames(cramer_matrix) <- colnames(mat_new_c8)
colnames(cramer_matrix) <- colnames(mat_new_c8) 
write.table(cramer_matrix, "top1920_all_cells_cramer_matrix.txt", quote = F, sep = "\t")

mat_new_c8_tmp <- mat_new_c8
mat_new_c8_tmp[mat_new_c8_tmp==2] <- 1
mat_new_c8_tmp[mat_new_c8_tmp==3] <- 1

cramer_matrix <- matrix(1:225,ncol = 15)
for (i in 1:15) {
  for (j in i:15) {
    cramer_matrix[i,j] <- cramerV(mat_new_c8_tmp[,i], mat_new_c8_tmp[,j])
    cramer_matrix[j,i] <- cramer_matrix[i,j]
  }
}
rownames(cramer_matrix) <- colnames(mat_new_c8_tmp)
colnames(cramer_matrix) <- colnames(mat_new_c8_tmp) 
write.table(cramer_matrix, "top1920_all_cells_cramer_matrix_only_active_repressive_un.txt", quote = F, sep = "\t")
