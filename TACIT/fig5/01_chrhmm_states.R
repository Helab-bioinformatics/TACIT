
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/")
files <- list.files("./")
files <- files[grep("bam", files)]
files <- as.data.frame(files)
library(tidyr)
files <- separate(files, files, c("v1", "v2", "v3", "v4", "v5"), sep = "_", remove = F)
files <- unite(files, ID, c("v1", "v2", "v3"), sep = "_")
write.table(files[,c("ID", "v4", "files")], "./cellmarkfiletable.txt", sep = "\t", quote = F, col.names = F, row.names = F)


setwd("./learnmodel_2kb/POSTERIOR/")
library(tidyr)
####add rownames for all posterior.txt
for (n in c(1:19, "X", "Y")) {
  post <- list.files("./")
  post <- post[grep(paste0("chr", n,"_"), post)]
  for (i in post) {
    dat <- read.csv(file = i , sep = "\t", header = F)
    dat <- dat[-1,]
    colnames(dat) <- dat[1,]
    dat <- dat[-1,]
    rownames(dat) <- paste0("chr",n,"_",1:nrow(dat))
    write.table(dat, file =paste0("./test/", i) , quote = F)
  }
}

####state-bin matrix for each cell
ID <- unique(files$ID)
for (i in ID) {
  data <- NULL
  for (j in c(1:19,"X", "Y")) {
    file <- paste0("./test/",i, "_12_chrhmm_chr", j, "_posterior.txt")
    data <- rbind(data, read.csv(file = file, header = T, sep = " "))
  }
  data$anno <- colnames(data)[max.col(data)]
  write.table(data, file=paste0("./test/", i, "_all_posterior.txt"), quote = F)
}


mat_1 <- read.csv(file = paste0("./test/",ID[1], "_all_posterior.txt"), sep = " ")
mat_2 <- read.csv(file = paste0("./test/",ID[2], "_all_posterior.txt"), sep = " ")
mat_3 <- read.csv(file = paste0("./test/",ID[3], "_all_posterior.txt"), sep = " ")
mat_4 <- read.csv(file = paste0("./test/",ID[4], "_all_posterior.txt"), sep = " ")
mat_5 <- read.csv(file = paste0("./test/",ID[5], "_all_posterior.txt"), sep = " ")
mat_6 <- read.csv(file = paste0("./test/",ID[6], "_all_posterior.txt"), sep = " ")
mat_7 <- read.csv(file = paste0("./test/",ID[7], "_all_posterior.txt"), sep = " ")
mat_8 <- read.csv(file = paste0("./test/",ID[8], "_all_posterior.txt"), sep = " ")
mat_9 <- read.csv(file = paste0("./test/",ID[9], "_all_posterior.txt"), sep = " ")
mat_10 <- read.csv(file = paste0("./test/",ID[10], "_all_posterior.txt"), sep = " ")
mat_11 <- read.csv(file = paste0("./test/",ID[11], "_all_posterior.txt"), sep = " ")
mat_12 <- read.csv(file = paste0("./test/",ID[12], "_all_posterior.txt"), sep = " ")
mat_13 <- read.csv(file = paste0("./test/",ID[13], "_all_posterior.txt"), sep = " ")

all_dat <- cbind(mat_1$anno, mat_2$anno, mat_3$anno, mat_4$anno, mat_5$anno, mat_6$anno, mat_7$anno, mat_8$anno, mat_9$anno, mat_10$anno,
                 mat_11$anno, mat_12$anno, mat_13$anno)
all_dat <- as.data.frame(all_dat)
rownames(all_dat) <- rownames(mat_1)
colnames(all_dat) <- ID
write.table(all_dat, "./test/all_cells_states.txt", quote = F)


##########################################################################################
library(tidyr)
library(matrixStats)
setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/learnmodel_2kb/POSTERIOR/")

all_mat <- read.csv("./test/all_cells_states.txt", sep = " ")
all_mat[all_mat=="E1"] <- 1
all_mat[all_mat=="E2"] <- 2
all_mat[all_mat=="E3"] <- 3
all_mat[all_mat=="E4"] <- 4
all_mat[all_mat=="E5"] <- 5
all_mat[all_mat=="E6"] <- 6
all_mat[all_mat=="E7"] <- 0
all_mat[all_mat=="E8"] <- 8
all_mat[all_mat=="E9"] <- 9
all_mat[all_mat=="E10"] <- 10
all_mat[all_mat=="E11"] <- 11
all_mat[all_mat=="E12"] <- 0
row <- rownames(all_mat)
all_mat <- apply(all_mat, 2, as.numeric)
rownames(all_mat) <- row
all_mat <- as.matrix(all_mat)
all_mat <- all_mat[which(rowSds(all_mat) >0),]
all_dat_2 <- all_mat[which(rowSums(all_mat == 0) <11),] # filtering

all_mat <- read.csv("./all_cells_states.txt", sep = " ")
all_dat_3 <- all_mat[rownames(all_dat_2),]
all_dat_3[all_dat_3=="E1"] <- 1
all_dat_3[all_dat_3=="E2"] <- 1
all_dat_3[all_dat_3=="E3"] <- 2
all_dat_3[all_dat_3=="E4"] <- 2
all_dat_3[all_dat_3=="E5"] <- 2
all_dat_3[all_dat_3=="E6"] <- 2
all_dat_3[all_dat_3=="E7"] <- 5
all_dat_3[all_dat_3=="E8"] <- 3
all_dat_3[all_dat_3=="E9"] <- 3
all_dat_3[all_dat_3=="E10"] <- 4
all_dat_3[all_dat_3=="E11"] <- 4
all_dat_3[all_dat_3=="E12"] <- 5

all_dat_3<- apply(all_dat_3, 2, as.numeric)
rownames(all_dat_3) <- rownames(all_dat_2)
all_dat_3 <- as.matrix(all_dat_3)
all_dat_3 <- all_dat_3[which(rowSds(all_dat_3) >0),]

library(pheatmap)
pdf("all_features_heatmap.pdf")
pheatmap(all_dat_3, 
         cluster_rows =T, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = F, 
         clustering_distance_rows="correlation",	
         clustering_method="ward.D2",
         color = colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

all_dat_3 <- all_dat_3[,c(1:11, 13)]
set.seed(15)
km<- kmeans(all_dat_3,5) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(all_dat_3, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column

library(pheatmap)
pdf("all_features_kmeans_heatmap.pdf")
pheatmap( m.kmeans[,1:ncol(m.kmeans)-1], 
          cluster_rows = F, cluster_cols = F, 
          show_rownames=FALSE, show_colnames=FALSE,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

m.kmeans <- as.data.frame(m.kmeans)
all_dat_3_r <- all_dat_3[rownames(m.kmeans[grep("4", m.kmeans$V13),]),]
all_dat_3_a <- all_dat_3[rownames(m.kmeans[grep("1|2|3|5", m.kmeans$V13),]),]
################################# randomForest ###############################
library(randomForest)
library(pROC)
library(ROCR)
library(matrixStats)
library(tidyr)
set.seed(12)

all_dat_3_select <- as.data.frame(t(all_dat_3_select))
all_dat_3_select$identity <- as.factor(c(rep("ICM", 6), rep("TE", 6)))
index <- sample(1:12, 3)#随机???3/10作为test，剩下为train
train <- all_dat_3_select[-index,]
test <- all_dat_3_select[index,]
err<-as.numeric()
set.seed(12)
#options(expressions = 5e5)
for(i in c(1:10)){
  mtry_test <- randomForest(identity ~., data=train,ntree=10000,mtry=i)
  err<- append( err, mean( mtry_test$err.rate ) )
}
mtry<-which.min(err)#mtry=6
set.seed(123)
mtry_fit <- randomForest(identity ~.,data = train, mtry=mtry, ntree=10000, importance=T, proximity=TRUE)
mtry_fit #mtry=9, ntree=10000, OOB error=11.11%

pdf("./test/train_pseudo_both_20230412_top780.pdf")
plot(mtry_fit)#best ntree=10000
plot(mtry_fit$err.rate)
title(main = "use all diff bins", sub = "err=33.33%")

train_predict <- predict(mtry_fit, train)
compare_train <- table(train_predict, train$identity)
compare_train
sum(diag(compare_train)/sum(compare_train))
test_predict <- predict(mtry_fit, test, type = "prob")
predicts <- t(apply(test_predict, 1, function(v){v/sum(v)}))
colnames(predicts) <- colnames(test_predict)
predicts <- data.frame(predicts, check.names = F)
predicts$predicted <- apply(predicts, 1, function(v){names(v)[max(v) == v]})
predicts$observed <- test$identity
#ROC
ROC <- pROC::roc(predicts$observed, as.numeric(predicts$ICM))
plot(ROC, print.auc=TRUE)

#get important features
importance_otu <- importance(mtry_fit)
varImpPlot(mtry_fit)
importance_otu <- as.data.frame(importance_otu)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)
#write.table(importance_otu, './randomForest/importance_otu_top1500.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(importance_otu[1:30, ], 'importance_otu_top30.txt', sep = '\t', col.names = NA, quote = FALSE)

# cross-validation
set.seed(123)
otu_train.cv <- replicate(5, rfcv(train[-ncol(train)], train$identity, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
#
library(ggplot2)
library(splines)  
pdf("test.pdf")
p <- ggplot(otu_train.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
p
dev.off()

###
all_dat_3_a <- all_dat_3_a[,rownames(importance_otu[1:600,])]
all_dat_3_r <- all_dat_3_r[,rownames(importance_otu[1:180,])]

##再跑一???

save(all_dat_3, all_dat_3_a,all_dat_3_r,mtry_fit, predict, all, cramer_matrix,test, train,importance_otu,file = "./chrhmm.random_20230412.Rdata")


setwd("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/04_chromHMM_pseudo/downsampling/learnmodel_2kb/POSTERIOR")
load("./test/chrhmm.random_20230412.Rdata")
################################ top features #################################
all_dat_3_select <- cbind(all_dat_3_a[,1:600], all_dat_3_r[,1:180])

mm10 <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/02_chromHMM_sc/learnmodel_2kb/POSTERIOR/test/randomForest/chrmm_loci.txt", sep = " ", header = F)
features <- mm10[which(mm10$V5 %in% colnames(all_dat_3_select)),1:4]
write.table(features[,2:4], "./test/20230412_top780.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# 
# ref <- read.csv("./test/gene.TEICM.tss50kb.region.bed", sep = "\t", header = F)
# ref$up <- ref$V2 -450000
# ref$down <- ref$V3 +450000
# write.table(ref[,c(1, 4:5)], "./gene.TEICM.tss500kb.region.bed", sep = "\t", quote = F, col.names = F, row.names = F)

######################## heatmap for blasto #################
library(pheatmap)
all_dat_3_select <- t(as.matrix(all_dat_3_select))
pdf("./test/all_features_heatmap_top780.pdf")
pheatmap(all_dat_3_select, 
         cluster_rows =T, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = F, 
         clustering_distance_rows="correlation",	
         clustering_method="ward.D2",
         color = colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()


set.seed(12)
km<- kmeans(all_dat_3_select,10) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(all_dat_3_select, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column

anno_row <- data.frame(cluster=m.kmeans[,ncol(m.kmeans)], cluster2=m.kmeans[,ncol(m.kmeans)])
rownames(anno_row) <- rownames(m.kmeans)

library(pheatmap)
pdf("./test/all_features_kmeans_top780_heatmap.pdf")
pheatmap( m.kmeans[,1:ncol(m.kmeans)-1], 
          cluster_rows = F, cluster_cols = F, 
          show_rownames=F, show_colnames=T,
          annotation_row = anno_row,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

m.kmeans <- as.data.frame(m.kmeans)
cluster_1 <- m.kmeans[which(m.kmeans$V13 == 1),]
cluster_2 <- m.kmeans[which(m.kmeans$V13 == 2),]
cluster_3 <- m.kmeans[which(m.kmeans$V13 == 3),]
cluster_4 <- m.kmeans[which(m.kmeans$V13 == 4),]
cluster_5 <- m.kmeans[which(m.kmeans$V13 == 5),]
cluster_6 <- m.kmeans[which(m.kmeans$V13 == 6),]
cluster_7 <- m.kmeans[which(m.kmeans$V13 == 7),]
cluster_8 <- m.kmeans[which(m.kmeans$V13 == 8),]
cluster_9 <- m.kmeans[which(m.kmeans$V13 == 9),]
cluster_10 <- m.kmeans[which(m.kmeans$V13 == 10),]
m.kmeans_new <- rbind(cluster_2, cluster_5, cluster_6, cluster_7, cluster_10,
                      cluster_1, cluster_3, cluster_4, cluster_8, cluster_9)
pdf("./test/all_features_kmeans_top780_heatmap_neworder.pdf")
pheatmap( m.kmeans_new[,1:ncol(m.kmeans_new)-1], 
          cluster_rows = F, cluster_cols = F, 
          show_rownames=F, show_colnames=T,
          annotation_row = anno_row,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

cluster_sample_1 <- rbind(cluster_2, cluster_5, cluster_6, cluster_7, cluster_10)
cluster_sample_1 <- cluster_sample_1[sample(nrow(cluster_sample_1)),]
cluster_sample_2 <- cluster_1[sample(nrow(cluster_1)),]
cluster_sample_3 <- rbind(cluster_3, cluster_4, cluster_8)
cluster_sample_3 <- cluster_sample_3[sample(nrow(cluster_sample_3)),]
cluster_sample_4 <- cluster_9[sample(nrow(cluster_9)),]
m.kmeans_new <- rbind(cluster_sample_1, cluster_sample_2, cluster_sample_3, cluster_sample_4)
pdf("./test/all_features_kmeans_top780_heatmap_neworder_sample.pdf")
pheatmap( m.kmeans_new[,1:ncol(m.kmeans_new)-1], 
          cluster_rows = F, cluster_cols = F, 
          show_rownames=F, show_colnames=T,
          annotation_row = anno_row,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

m.kmeans_new <- rbind(cluster_sample_1, cluster_1, cluster_sample_3, cluster_9)
pdf("./test/all_features_kmeans_top780_heatmap_neworder_sample_2.pdf")
pheatmap( m.kmeans_new[,1:ncol(m.kmeans_new)-1], 
          cluster_rows = F, cluster_cols = F, 
          show_rownames=F, show_colnames=T,
          annotation_row = anno_row,
          color =colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

######## for TE-specific regions #########
active <- mm10[which(mm10$V5 %in% rownames(cluster_sample_1)),]
active_2 <- mm10[which(mm10$V5 %in% rownames(cluster_1[1:8,])),]
active <- rbind(active, active_2)
repressive <- mm10[which(mm10$V5 %in% rownames(cluster_1[8:nrow(cluster_1),])),]
write.table(active[,2:4], "./TE_specific_active_chromatin_regions.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(repressive[,2:4], "./TE_specific_repressive_chromatin_regions.bed", sep = "\t", quote = F, col.names = F, row.names = F)

TE_specific <- rbind(cluster_sample_1, cluster_1)
a <- as.data.frame(table(TE_specific$TE_pseudo_7))

TE <- read.csv("./test/TE_specific_ratio.csv", sep = ",")
TE <- TE[grep("5", invert = T, TE$type),]
pdf("./test/TE_specific_region_ratio.pdf")
ggplot(TE, aes(x=cell, y=fre))+
  geom_smooth(aes(group=type), method = "loess", se = FALSE)+
  theme_bw()
dev.off()


######## for ICM-specific regions #########
active <- mm10[which(mm10$V5 %in% rownames(cluster_sample_3)),]
active_2 <- mm10[which(mm10$V5 %in% rownames(cluster_9[1:49,])),]
active <- rbind(active, active_2)
repressive <- mm10[which(mm10$V5 %in% rownames(cluster_9[50:nrow(cluster_1),])),]
write.table(active[,2:4], "./ICM_specific_active_chromatin_regions.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(repressive[,2:4], "./ICM_specific_repressive_chromatin_regions.bed", sep = "\t", quote = F, col.names = F, row.names = F)

ICM_specific <- rbind(cluster_sample_3, cluster_9)
a <- as.data.frame(table(ICM_specific$ICM_pseudo_1))

ICM <- read.csv("./test/ICM_specific_ratio.csv", sep = ",")
ICM <- ICM[grep("5", invert = T, ICM$type),]
pdf("./test/ICM_specific_region_ratio.pdf")
ggplot(ICM, aes(x=cell, y=fre))+
  geom_smooth(aes(group=type), method = "loess", se = FALSE)+
  theme_bw()
dev.off()


dat <- data.frame(group=c("proximal", "distal"), value=c(486, 294))
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
pdf("./test/top780_around_ICMTE_marker_genes_500kb.pdf.pdf")
ggplot(dat, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start=0)+
  scale_fill_brewer("Blues") + blank_theme +
  theme(axis.text.x=element_blank())
dev.off()


########################## predict #######################
library(randomForest)
library(pROC)
library(ROCR)
library(matrixStats)
library(tidyr)
set.seed(12)

pre_dat <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/06_blasto_inte/03_scChromHMM/05_predict_pseudo/downsampling/learnmodel_2kb/POSTERIOR/test/all_cells_states.txt", sep = " ")
pre_dat <- pre_dat[rownames(all_dat_3_select),]
pre_dat[pre_dat=="E1"] <- 1
pre_dat[pre_dat=="E2"] <- 1
pre_dat[pre_dat=="E3"] <- 2
pre_dat[pre_dat=="E4"] <- 2
pre_dat[pre_dat=="E5"] <- 2
pre_dat[pre_dat=="E6"] <- 2
pre_dat[pre_dat=="E7"] <- 5
pre_dat[pre_dat=="E8"] <- 3
pre_dat[pre_dat=="E9"] <- 3
pre_dat[pre_dat=="E10"] <- 4
pre_dat[pre_dat=="E11"] <- 4
pre_dat[pre_dat=="E12"] <- 5

pre_dat <- t(pre_dat)
predict <- predict(mtry_fit,pre_dat, type = "prob")
predict <- as.data.frame(predict)

predict[which(predict$ICM>0.6),3] <- "ICM"
predict[which(predict$TE>0.6),3] <- "TE"

morula <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/20230105_allStage_integration/05_morula_inte/based_cluster/downsampling_4/downsampling_3/learnmodel_2kb/POSTERIOR/test/all_cells_states.txt", sep = " ")
morula <- morula[rownames(all_dat_3_select),]
morula[morula=="E1"] <- 1
morula[morula=="E2"] <- 1
morula[morula=="E3"] <- 2
morula[morula=="E4"] <- 2
morula[morula=="E5"] <- 2
morula[morula=="E6"] <- 2
morula[morula=="E7"] <- 5
morula[morula=="E8"] <- 3
morula[morula=="E9"] <- 3
morula[morula=="E10"] <- 4
morula[morula=="E11"] <- 4
morula[morula=="E12"] <- 5

morula <- t(morula)
predict_morula <- predict(mtry_fit,morula, type = "prob")
predict_morula <- as.data.frame(predict_morula)

morula <- morula[-3,]
######################## similarity matrix #####################
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

blasto <- as.data.frame(t(all_dat_3_select))
blasto <- blasto[!(rownames(blasto) %in% c("TE_pseudo_2", "ICM_pseudo_5")),]
pre_dat <- as.data.frame(pre_dat)
pre_dat <- pre_dat[1:12,]
all <- rbind(pre_dat, blasto, morula)
all <- as.data.frame(t(all))
all <- all[,c("X4cell_2", "X8cell_1", "X8cell_2", 
              "X4cell_3", "X4cell_4", "X8cell_5", "X8cell_6", "X8cell_7", "morula_c10c12","morula_c1c15","morula_c4c6","morula_c4c9", "ICM_pseudo_1", "ICM_pseudo_3", "ICM_pseudo_2", "ICM_pseudo_4", "ICM_pseudo_6",
              "X4cell_1", "X8cell_3", "X8cell_4", "X8cell_8", "morula_c5c7", "TE_pseudo_1", "TE_pseudo_3", "TE_pseudo_4", "TE_pseudo_5", "TE_pseudo_7")]
# 计算所有列之间的CramerV相关???

cramer_matrix <- matrix(1:729,ncol = 27)
for (i in 1:27) {
  for (j in i:27) {
    cramer_matrix[i,j] <- cramerV(all[,i], all[,j])
    cramer_matrix[j,i] <- cramer_matrix[i,j]
  }
}
rownames(cramer_matrix) <- colnames(all)
colnames(cramer_matrix) <- colnames(all)

pdf("./test/all_features_heatmap_top780_4cellToblasto_cramerV_20230420.pdf")
cramer_matrix[cramer_matrix==1] <- 0
pheatmap(cramer_matrix, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = T, scale ="none",
         show_colnames = T, 
         colorRampPalette(c("#2a628c","#3d8ec1","#f2ee75","#f9fa02"))(500))
dev.off()

all <- apply(all, 2, as.numeric)
rownames(all) <- rownames(all_dat_3_select)
cor_pearson <- cor(all, method = 'pearson')
cor_spearman <- cor(all, method = 'spearman')

pdf("./test/all_features_heatmap_top780_4cellToblasto_spearman_20230420.pdf")
cor_spearman[cor_spearman==1] <- 0
pheatmap(cor_spearman, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = T, scale ="none",
         show_colnames = T, 
         colorRampPalette(c("#2a628c","#3d8ec1","#f2ee75","#f9fa02"))(500))
dev.off()


all <- all[rownames(m.kmeans_new),]
pdf("./test/all_features_heatmap_top780_4cellToblasto_heatmap_20230420.pdf")
pheatmap(all, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = T, 
         color = colorRampPalette(c( "#117733", "#f4da45", "#db0b20", "#169dd3","#fafafc" ))(5))
dev.off()

################################# ICM and TE specific ###############################
all_ICM <- all[rownames(ICM_specific),]
all_ICM_1 <- all_ICM[, c("X4cell_3", "X4cell_4", "X8cell_5", "X8cell_6", "X8cell_7", "morula_c10c12","morula_c1c15","morula_c4c6","morula_c4c9", "ICM_pseudo_1", "ICM_pseudo_3", "ICM_pseudo_2", "ICM_pseudo_4", "ICM_pseudo_6")]
fre_ICM_1 <- as.data.frame(table(all_ICM_1[,1]))
for (i in 2:14) {
  a <- as.data.frame(table(all_ICM_1[,i]))
  fre_ICM_1[,i+1] <- a[,2]
}
colnames(fre_ICM_1) <- c("state", "X4cell_3", "X4cell_4", "X8cell_5", "X8cell_6", "X8cell_7", "morula_c10c12","morula_c1c15","morula_c4c6","morula_c4c9", "ICM_pseudo_1", "ICM_pseudo_3", "ICM_pseudo_2", "ICM_pseudo_4", "ICM_pseudo_6")
fre_ICM_1 <- gather(fre_ICM_1, key="cells",value="fre", -state)
fre_ICM_1 <- fre_ICM_1[grep("5", fre_ICM_1$state, invert = T),]
fre_ICM_1$cells <- factor(fre_ICM_1$cells, levels = c("X4cell_3", "X4cell_4", "X8cell_5", "X8cell_6", "X8cell_7", "morula_c10c12","morula_c1c15","morula_c4c6","morula_c4c9", "ICM_pseudo_1", "ICM_pseudo_3", "ICM_pseudo_2", "ICM_pseudo_4", "ICM_pseudo_6"))

pdf("./test/ICM_lineage_state_ratio.pdf")
ggplot( fre_ICM_1, aes( x = cells, weight = fre, fill = state))+
  geom_bar( position = "stack")+
  theme_bw()

ggplot(fre_ICM_1, aes(x=cells, y=fre))+
  geom_smooth(aes(group=state), method = "loess", se = FALSE)+
  theme_bw()
dev.off()
write.table(fre_ICM_1, "./test/ICM_lineage_state_ratio.txt", sep = "\t", quote = F)


all_TE <- all[rownames(TE_specific),]
all_TE_1 <- all_TE[, c("X4cell_1", "X8cell_3", "X8cell_4", "X8cell_8", "morula_c5c7", "TE_pseudo_1", "TE_pseudo_3", "TE_pseudo_4", "TE_pseudo_5", "TE_pseudo_7")]

fre_TE_1 <- as.data.frame(table(all_TE_1[,1]))
for (i in 2:10) {
  a <- as.data.frame(table(all_TE_1[,i]))
  fre_TE_1[,i+1] <- a[,2]
}
colnames(fre_TE_1) <- c("state", "X4cell_1", "X8cell_3", "X8cell_4", "X8cell_8", "morula_c5c7", "TE_pseudo_1", "TE_pseudo_3", "TE_pseudo_4", "TE_pseudo_5", "TE_pseudo_7")
fre_TE_1 <- gather(fre_TE_1, key="cells",value="fre", -state)
fre_TE_1 <- fre_TE_1[grep("5", fre_TE_1$state, invert = T),]
fre_TE_1$cells <- factor(fre_TE_1$cells, levels = c("X4cell_1", "X8cell_3", "X8cell_4", "X8cell_8", "morula_c5c7", "TE_pseudo_1", "TE_pseudo_3", "TE_pseudo_4", "TE_pseudo_5", "TE_pseudo_7"))

pdf("./test/TE_lineage_state_ratio.pdf")
ggplot( fre_TE_1, aes( x = cells, weight = fre, fill = state))+
  geom_bar( position = "stack")+
  theme_bw()

ggplot(fre_TE_1, aes(x=cells, y=fre))+
  geom_smooth(aes(group=state), method = "loess", se = FALSE)+
  theme_bw()
dev.off()
write.table(fre_TE_1, "./test/TE_lineage_state_ratio.txt", sep = "\t", quote = F)
save(all_dat_3, all_dat_3_a, all_dat_3_r, all_dat_3_select, cramer_matrix, m.kmeans, m.kmeans_new, mm10,mtry_fit, predict, all, ICM_specific, all_ICM,all_TE,TE_specific,file = "./test/20230413_chrhmm.Rdata")


########################## differential regions for calling TF motifs ####################
########### for TE ###########
all_TE <- as.data.frame(all_TE)
TE_4cell <- all_TE[,grep("4cell_1", colnames(all_TE))]
TE_4cell_en <- rownames(all_TE)[TE_4cell ==2]
TE_4cell_pr <- rownames(all_TE)[TE_4cell ==1]
TE_4cell_tr <- rownames(all_TE)[TE_4cell ==3]
TE_4cell_re <- rownames(all_TE)[TE_4cell ==4]

TE_4cell_en <- mm10[which(mm10$V5 %in% TE_4cell_en),]
TE_4cell_pr <- mm10[which(mm10$V5 %in% TE_4cell_pr),]
TE_4cell_tr <- mm10[which(mm10$V5 %in% TE_4cell_tr),]
TE_4cell_re <- mm10[which(mm10$V5 %in% TE_4cell_re),]
TE_4cell_active <- rbind(TE_4cell_en, TE_4cell_pr, TE_4cell_tr)

write.table(TE_4cell_en[,c("V2", "V3", "V4")], "./test/states/TE_4cell_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_4cell_re[,c("V2", "V3", "V4")], "./test/states/TE_4cell_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_4cell_active[,c("V2", "V3", "V4")], "./test/states/TE_4cell_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)



TE_8cell <- all_TE[,grep("8cell_3|8cell_4|8cell_8", colnames(all_TE))]
TE_8cell_en <- rownames(all_TE)[TE_8cell ==2]
TE_8cell_pr <- rownames(all_TE)[TE_8cell ==1]
TE_8cell_tr <- rownames(all_TE)[TE_8cell ==3]
TE_8cell_re <- rownames(all_TE)[TE_8cell ==4]

TE_8cell_en <- mm10[which(mm10$V5 %in% TE_8cell_en),]
TE_8cell_pr <- mm10[which(mm10$V5 %in% TE_8cell_pr),]
TE_8cell_tr <- mm10[which(mm10$V5 %in% TE_8cell_tr),]
TE_8cell_re <- mm10[which(mm10$V5 %in% TE_8cell_re),]
TE_8cell_active <- rbind(TE_8cell_en, TE_8cell_pr, TE_8cell_tr)

write.table(TE_8cell_en[,c("V2", "V3", "V4")], "./test/states/TE_8cell_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_8cell_re[,c("V2", "V3", "V4")], "./test/states/TE_8cell_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_8cell_active[,c("V2", "V3", "V4")], "./test/states/TE_8cell_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)



TE_morula <- all_TE[,grep("morula_c5c7", colnames(all_TE))]
TE_morula_en <- rownames(all_TE)[TE_morula ==2]
TE_morula_pr <- rownames(all_TE)[TE_morula ==1]
TE_morula_tr <- rownames(all_TE)[TE_morula ==3]
TE_morula_re <- rownames(all_TE)[TE_morula ==4]

TE_morula_en <- mm10[which(mm10$V5 %in% TE_morula_en),]
TE_morula_pr <- mm10[which(mm10$V5 %in% TE_morula_pr),]
TE_morula_tr <- mm10[which(mm10$V5 %in% TE_morula_tr),]
TE_morula_re <- mm10[which(mm10$V5 %in% TE_morula_re),]
TE_morula_active <- rbind(TE_morula_en, TE_morula_pr, TE_morula_tr)

write.table(TE_morula_en[,c("V2", "V3", "V4")], "./test/states/TE_morula_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_morula_re[,c("V2", "V3", "V4")], "./test/states/TE_morula_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_morula_active[,c("V2", "V3", "V4")], "./test/states/TE_morula_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)



TE_blasto <- all_TE[,grep("TE", colnames(all_TE))]
TE_blasto_en <- rownames(all_TE)[TE_blasto ==2]
TE_blasto_pr <- rownames(all_TE)[TE_blasto ==1]
TE_blasto_tr <- rownames(all_TE)[TE_blasto ==3]
TE_blasto_re <- rownames(all_TE)[TE_blasto ==4]

TE_blasto_en <- mm10[which(mm10$V5 %in% TE_blasto_en),]
TE_blasto_pr <- mm10[which(mm10$V5 %in% TE_blasto_pr),]
TE_blasto_tr <- mm10[which(mm10$V5 %in% TE_blasto_tr),]
TE_blasto_re <- mm10[which(mm10$V5 %in% TE_blasto_re),]
TE_blasto_active <- rbind(TE_blasto_en, TE_blasto_pr, TE_blasto_tr)

write.table(TE_blasto_en[,c("V2", "V3", "V4")], "./test/states/TE_blasto_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_blasto_re[,c("V2", "V3", "V4")], "./test/states/TE_blasto_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(TE_blasto_active[,c("V2", "V3", "V4")], "./test/states/TE_blasto_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)


########### for ICM ###########
all_ICM <- as.data.frame(all_ICM)
ICM_4cell <- all_ICM[,grep("4cell_3|4cell_4", colnames(all_ICM))]
ICM_4cell_en <- rownames(all_ICM)[ICM_4cell ==2]
ICM_4cell_pr <- rownames(all_ICM)[ICM_4cell ==1]
ICM_4cell_tr <- rownames(all_ICM)[ICM_4cell ==3]
ICM_4cell_re <- rownames(all_ICM)[ICM_4cell ==4]

ICM_4cell_en <- mm10[which(mm10$V5 %in% ICM_4cell_en),]
ICM_4cell_pr <- mm10[which(mm10$V5 %in% ICM_4cell_pr),]
ICM_4cell_tr <- mm10[which(mm10$V5 %in% ICM_4cell_tr),]
ICM_4cell_re <- mm10[which(mm10$V5 %in% ICM_4cell_re),]
ICM_4cell_active <- rbind(ICM_4cell_en, ICM_4cell_pr, ICM_4cell_tr)

write.table(ICM_4cell_en[,c("V2", "V3", "V4")], "./test/states/ICM_4cell_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_4cell_re[,c("V2", "V3", "V4")], "./test/states/ICM_4cell_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_4cell_active[,c("V2", "V3", "V4")], "./test/states/ICM_4cell_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)



ICM_8cell <- all_ICM[,grep("8cell_5|8cell_6|8cell_7", colnames(all_ICM))]
ICM_8cell_en <- rownames(all_ICM)[ICM_8cell ==2]
ICM_8cell_pr <- rownames(all_ICM)[ICM_8cell ==1]
ICM_8cell_tr <- rownames(all_ICM)[ICM_8cell ==3]
ICM_8cell_re <- rownames(all_ICM)[ICM_8cell ==4]

ICM_8cell_en <- mm10[which(mm10$V5 %in% ICM_8cell_en),]
ICM_8cell_pr <- mm10[which(mm10$V5 %in% ICM_8cell_pr),]
ICM_8cell_tr <- mm10[which(mm10$V5 %in% ICM_8cell_tr),]
ICM_8cell_re <- mm10[which(mm10$V5 %in% ICM_8cell_re),]
ICM_8cell_active <- rbind(ICM_8cell_en, ICM_8cell_pr, ICM_8cell_tr)

write.table(ICM_8cell_en[,c("V2", "V3", "V4")], "./test/states/ICM_8cell_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_8cell_re[,c("V2", "V3", "V4")], "./test/states/ICM_8cell_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_8cell_active[,c("V2", "V3", "V4")], "./test/states/ICM_8cell_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)



ICM_morula <- all_ICM[,grep("c10c12|c1c15|c4c6|c4c9", colnames(all_ICM))]
ICM_morula_en <- rownames(all_ICM)[ICM_morula ==2]
ICM_morula_pr <- rownames(all_ICM)[ICM_morula ==1]
ICM_morula_tr <- rownames(all_ICM)[ICM_morula ==3]
ICM_morula_re <- rownames(all_ICM)[ICM_morula ==4]

ICM_morula_en <- mm10[which(mm10$V5 %in% ICM_morula_en),]
ICM_morula_pr <- mm10[which(mm10$V5 %in% ICM_morula_pr),]
ICM_morula_tr <- mm10[which(mm10$V5 %in% ICM_morula_tr),]
ICM_morula_re <- mm10[which(mm10$V5 %in% ICM_morula_re),]
ICM_morula_active <- rbind(ICM_morula_en, ICM_morula_pr, ICM_morula_tr)

write.table(ICM_morula_en[,c("V2", "V3", "V4")], "./test/states/ICM_morula_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_morula_re[,c("V2", "V3", "V4")], "./test/states/ICM_morula_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_morula_active[,c("V2", "V3", "V4")], "./test/states/ICM_morula_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)



ICM_blasto <- all_ICM[,grep("ICM", colnames(all_ICM))]
ICM_blasto_en <- rownames(all_ICM)[ICM_blasto ==2]
ICM_blasto_pr <- rownames(all_ICM)[ICM_blasto ==1]
ICM_blasto_tr <- rownames(all_ICM)[ICM_blasto ==3]
ICM_blasto_re <- rownames(all_ICM)[ICM_blasto ==4]

ICM_blasto_en <- mm10[which(mm10$V5 %in% ICM_blasto_en),]
ICM_blasto_pr <- mm10[which(mm10$V5 %in% ICM_blasto_pr),]
ICM_blasto_tr <- mm10[which(mm10$V5 %in% ICM_blasto_tr),]
ICM_blasto_re <- mm10[which(mm10$V5 %in% ICM_blasto_re),]
ICM_blasto_active <- rbind(ICM_blasto_en, ICM_blasto_pr, ICM_blasto_tr)

write.table(ICM_blasto_en[,c("V2", "V3", "V4")], "./test/states/ICM_blasto_enhancer.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_blasto_re[,c("V2", "V3", "V4")], "./test/states/ICM_blasto_repressive.bed", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(ICM_blasto_active[,c("V2", "V3", "V4")], "./test/states/ICM_blasto_active.bed", sep = "\t", quote = F, col.names = F, row.names = F)

