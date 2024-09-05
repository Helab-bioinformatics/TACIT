
##################################################################################################
###################################### kmeans for coclustering ###################################
setwd("/media/helab/data3/min/TACIT/coTarget_early2cell/cicero_h3k4me3_h3k27ac/all_tss5kb/ZGA_related_cutoff0.2")
up <- read.csv("./up_promoter_enhancer_h3k4me3.txt", sep = "\t")
down <- read.csv("./down_promoter_enhancer_h3k4me3.txt", sep = "\t")
shared <- read.csv("./shared_promoter_enhancer_h3k4me3.txt", sep = "\t")
up_k27 <- read.csv("./up_promoter_enhancer_h3k27ac.txt", sep = "\t")
down_k27 <- read.csv("./down_promoter_enhancer_h3k27ac.txt", sep = "\t")
shared_k27 <- read.csv("./shared_promoter_enhancer_h3k27ac.txt", sep = "\t")


### depth for normalization
k4me3_id <- read.csv("/media/helab/data3/min/TACIT/coTarget_early2cell/01_h3k4me3/rename/rmdupID.txt", header = F)
k4me3_depth <- read.csv("/media/helab/data3/min/TACIT/coTarget_early2cell/01_h3k4me3/rename/rmdup.txt", header = F)
k27ac_id <- read.csv("/media/helab/data3/min/TACIT/coTarget_early2cell/02_h3k27ac/rename/rmdupID.txt", header = F)
k27ac_depth <- read.csv("/media/helab/data3/min/TACIT/coTarget_early2cell/02_h3k27ac/rename/rmdup.txt", header = F)

k4me3_depth <- data.frame(cells=k4me3_id, depth=k4me3_depth)
k27ac_depth <- data.frame(cells=k27ac_id, depth=k27ac_depth)
rownames(k4me3_depth) <- paste0(k4me3_depth$V1, ".bam")
rownames(k27ac_depth) <- paste0(k27ac_depth$V1, ".bam")

### normalize peak length
library(tidyr)
up <- separate(up, loci, c("loci", "min"), sep = ":")
up <- separate(up, min, c("up", "down"), sep = "-")
up$length <- as.numeric(up$down) - as.numeric(up$up)

up_k27 <- separate(up_k27, loci, c("loci", "min"), sep = ":")
up_k27 <- separate(up_k27, min, c("up", "down"), sep = "-")
up_k27$length <- as.numeric(up_k27$down) - as.numeric(up_k27$up)

down <- separate(down, loci, c("loci", "min"), sep = ":")
down <- separate(down, min, c("up", "down"), sep = "-")
down$length <- as.numeric(down$down) - as.numeric(down$up)

down_k27 <- separate(down_k27, loci, c("loci", "min"), sep = ":")
down_k27 <- separate(down_k27, min, c("up", "down"), sep = "-")
down_k27$length <- as.numeric(down_k27$down) - as.numeric(down_k27$up)

shared <- separate(shared, loci, c("loci", "min"), sep = ":")
shared <- separate(shared, min, c("up", "down"), sep = "-")
shared$length <- as.numeric(shared$down) - as.numeric(shared$up)

shared_k27 <- separate(shared_k27, loci, c("loci", "min"), sep = ":")
shared_k27 <- separate(shared_k27, min, c("up", "down"), sep = "-")
shared_k27$length <- as.numeric(shared_k27$down) - as.numeric(shared_k27$up)

#########################################################################
##################################### up #################################
#### h3k4me3 ####
up_k4 <- up[,4:93]
up_k4 <- apply(up_k4, 2, as.numeric)
rownames(up_k4) <- rownames(up)

## normalize depth
k4me3_depth <- k4me3_depth[colnames(up_k4),]
up_k4_norde <- sweep(up_k4, 2, k4me3_depth$V1.1, "/")
## normalize peak length
up_k4_nor <- sweep(up_k4_norde, 1, up$length, "/")

dat1 <- matrix(1:78522,ncol=69)#1138
n=1
for (i in 1:69) {
  m=i+20
  dat1[,n] <- rowSums(up_k4_nor[,i:m])
  n=n+1
}
rownames(dat1) <- rownames(up)

#### h3k27ac ####
up_reads <- up_k27[,4:92]
up_reads <- apply(up_reads, 2, as.numeric)
rownames(up_reads) <- rownames(up_k27)

## normalize depth
k27ac_depth <- k27ac_depth[colnames(up_reads),]
up_k27_norde <- sweep(up_reads, 2, k27ac_depth$V1.1, "/")
## normalize peak length
up_k27_nor <- sweep(up_k27_norde, 1, up_k27$length, "/")

dat2 <- matrix(1:78522,ncol=69)#1138
n=1
for (i in 1:69) {
  m=i+20
  dat2[,n] <- rowSums(up_k27_nor[,i:m])
  n=n+1
}
rownames(dat2) <- rownames(up_k27_nor)

######## Z-score #####
mean_value <- mean(dat1)
std_dev <- sd(as.vector(dat1))
z_score_matrix_1 <- (dat1 - mean_value) / std_dev

mean_value <- mean(dat2)
std_dev <- sd(as.vector(dat2))
z_score_matrix_2 <- (dat2 - mean_value) / std_dev

dat1_tmp2 <- scale(dat1)
dat1_tmp2[dat1_tmp2>1.5] <- 1.5
dat2_tmp2 <- scale(dat2)
dat2_tmp2[dat2_tmp2>1.5] <- 1.5

dat <- cbind(z_score_matrix_1, z_score_matrix_2)
dat[dat>1.5] <- 1.5

set.seed(11)
km<- kmeans(dat,2) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(dat, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column

pdf("./ZGA_zscore_kmeans/up_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_zscore_normalize.pdf")
pheatmap( m.kmeans[,1:ncol(m.kmeans)-1], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

z_score_matrix_1 <- z_score_matrix_1[rownames(m.kmeans),]
z_score_matrix_1[z_score_matrix_1 >2.5] <- 2.5
pdf("./ZGA_zscore_kmeans/up_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore_normalize.pdf")
pheatmap( z_score_matrix_1, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 

z_score_matrix_2 <- z_score_matrix_2[rownames(m.kmeans),]
z_score_matrix_2[z_score_matrix_2 >2] <- 2
pdf("./ZGA_zscore_kmeans/up_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore_normalize.pdf")
pheatmap( z_score_matrix_2, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

############### curve ###################
z_score_matrix_1_c1 <- z_score_matrix_1[1:649,]
z_score_matrix_1_c2 <- z_score_matrix_1[650:1136,]
z_score_matrix_2_c1 <- z_score_matrix_2[1:649,]
z_score_matrix_2_c2 <- z_score_matrix_2[650:1136,]

dat <- data.frame(k4me3_c1=apply(z_score_matrix_1_c1, 2, median), k4me3_c2=apply(z_score_matrix_1_c2, 2, median),
                  k27ac_c1=apply(z_score_matrix_2_c1, 2, median), k27ac_c2=apply(z_score_matrix_2_c2, 2, median),
                  pseuo=1:69)
library(ggplot2)
pdf("./ZGA_zscore_kmeans/up_promoter_enhancer_cocluster_curve_k4k27_zscore.pdf")
ggplot(dat, aes(x=pseuo))+
  geom_smooth(aes(y=k4me3_c1), method = "loess", color = "#9bd3b1", se = FALSE)+
  geom_smooth(aes(y=k4me3_c2), method = "loess", color = "#17964d", se = FALSE)+
  geom_smooth(aes(y=k27ac_c1), method = "loess", color = "#f4af9d", se = FALSE)+
  geom_smooth(aes(y=k27ac_c2), method = "loess", color = "#f15a29", se = FALSE)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #????x???̶ȱ?ǩ????????ʾ??б?Ƕ?Ϊ15?ȣ??????µ???1(hjust = 1)????????ΪTimes??СΪ20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #????y???̶ȱ?ǩ???????أ???????С????????ʽΪplain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #????y????????????????
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #ȥ??Ĭ???????Ļ?ɫ??????x=0????y=0???Ӵ???ʾ(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #????ͼ?????ӱ?????????????
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #????ͼ?????ܱ?????????????
                                  size=12),
        panel.grid.major = element_blank())  #????ʾ??????
dev.off()

############################################################################
##################################### down #################################
#### h3k4me3 ####
down_k4 <- down[,4:93]
down_k4 <- apply(down_k4, 2, as.numeric)
rownames(down_k4) <- rownames(down)

## normalize depth
k4me3_depth <- k4me3_depth[colnames(down_k4),]
down_k4_norde <- sweep(down_k4, 2, k4me3_depth$V1.1, "/")
## normalize peak length
down_k4_nor <- sweep(down_k4_norde, 1, down$length, "/")

dat1 <- matrix(1:35535,ncol=69)#515
n=1
for (i in 1:69) {
  m=i+20
  dat1[,n] <- rowSums(down_k4_nor[,i:m])
  n=n+1
}
rownames(dat1) <- rownames(down)

#### h3k27ac ####
down_reads <- down_k27[,4:92]
down_reads <- apply(down_reads, 2, as.numeric)
rownames(down_reads) <- rownames(down_k27)

## normalize depth
k27ac_depth <- k27ac_depth[colnames(down_reads),]
down_k27_norde <- sweep(down_reads, 2, k27ac_depth$V1.1, "/")
## normalize peak length
down_k27_nor <- sweep(down_k27_norde, 1, down_k27$length, "/")

dat2 <- matrix(1:35535,ncol=69)#515
n=1
for (i in 1:69) {
  m=i+20
  dat2[,n] <- rowSums(down_k27_nor[,i:m])
  n=n+1
}
rownames(dat2) <- rownames(down_k27_nor)

######## z-score #####
mean_value <- mean(dat1)
std_dev <- sd(as.vector(dat1))
z_score_matrix_1 <- (dat1 - mean_value) / std_dev

mean_value <- mean(dat2)
std_dev <- sd(as.vector(dat2))
z_score_matrix_2 <- (dat2 - mean_value) / std_dev

dat1_tmp2 <- scale(dat1)
dat1_tmp2[dat1_tmp2>1.5] <- 1.5
dat2_tmp2 <- scale(dat2)
dat2_tmp2[dat2_tmp2>1.5] <- 1.5

dat <- cbind(z_score_matrix_1, z_score_matrix_2)
dat[dat>1.5] <- 1.5

set.seed(11)
km<- kmeans(dat,2) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(dat, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column

pdf("./ZGA_zscore_kmeans/down_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_zscore_normalize.pdf")
pheatmap( m.kmeans[,1:ncol(m.kmeans)-1], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

z_score_matrix_1 <- z_score_matrix_1[rownames(m.kmeans),]
z_score_matrix_1[z_score_matrix_1 >2.5] <- 2.5
pdf("./ZGA_zscore_kmeans/down_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore_normalize.pdf")
pheatmap( z_score_matrix_1, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 

z_score_matrix_2 <- z_score_matrix_2[rownames(m.kmeans),]
z_score_matrix_2[z_score_matrix_2 >2] <- 2
pdf("./ZGA_zscore_kmeans/down_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore_normalize.pdf")
pheatmap( z_score_matrix_2, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

############### curve ###################
z_score_matrix_1_c1 <- z_score_matrix_1[1:182,]
z_score_matrix_1_c2 <- z_score_matrix_1[183:513,]
z_score_matrix_2_c1 <- z_score_matrix_2[1:182,]
z_score_matrix_2_c2 <- z_score_matrix_2[183:513,]

dat <- data.frame(k4me3_c1=apply(z_score_matrix_1_c1, 2, median), k4me3_c2=apply(z_score_matrix_1_c2, 2, median),
                  k27ac_c1=apply(z_score_matrix_2_c1, 2, median), k27ac_c2=apply(z_score_matrix_2_c2, 2, median),
                  pseuo=1:69)
library(ggplot2)
pdf("./ZGA_zscore_kmeans/down_promoter_enhancer_cocluster_curve_k4k27_zscore.pdf")
ggplot(dat, aes(x=pseuo))+
  geom_smooth(aes(y=k4me3_c1), method = "loess", color = "#9bd3b1", se = FALSE)+
  geom_smooth(aes(y=k4me3_c2), method = "loess", color = "#17964d", se = FALSE)+
  geom_smooth(aes(y=k27ac_c1), method = "loess", color = "#f4af9d", se = FALSE)+
  geom_smooth(aes(y=k27ac_c2), method = "loess", color = "#f15a29", se = FALSE)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #????x???̶ȱ?ǩ????????ʾ??б?Ƕ?Ϊ15?ȣ??????µ???1(hjust = 1)????????ΪTimes??СΪ20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #????y???̶ȱ?ǩ???????أ???????С????????ʽΪplain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #????y????????????????
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #ȥ??Ĭ???????Ļ?ɫ??????x=0????y=0???Ӵ???ʾ(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #????ͼ?????ӱ?????????????
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #????ͼ?????ܱ?????????????
                                  size=12),
        panel.grid.major = element_blank())  #????ʾ??????
dev.off()

#################################### shared #################################
shared <- read.csv("./shared_promoter_enhancer_h3k4me3.txt", sep = "\t")
shared_k27 <- read.csv("./shared_promoter_enhancer_h3k27ac.txt", sep = "\t")

shared_k4 <- shared[,2:90]
shared_k4 <- apply(shared_k4, 2, as.numeric)
rownames(shared_k4) <- rownames(shared)
dat1 <- matrix(1:10971,ncol=69)#1159
n=1
for (i in 1:69) {
  m=i+20
  dat1[,n] <- rowSums(shared_k4[,i:m])
  n=n+1
}
rownames(dat1) <- rownames(shared)

shared_reads <- shared_k27[,2:90]
shared_reads <- apply(shared_reads, 2, as.numeric)
rownames(shared_reads) <- rownames(shared_k27)
dat2 <- matrix(1:10971,ncol=69)#159
n=1
for (i in 1:69) {
  m=i+20
  dat2[,n] <- rowSums(shared_reads[,i:m])
  n=n+1
}
rownames(dat2) <- rownames(shared_k27)

dat1_tmp <- t(dat1)
dat1_tmp <- scale(dat1_tmp)
dat1_tmp <- t(dat1_tmp)
dat1_tmp[dat1_tmp>4] <- 4
dat2_tmp <- t(dat2)
dat2_tmp <- scale(dat2_tmp)
dat2_tmp <- t(dat2_tmp)
dat2_tmp[dat2_tmp>4] <- 4

dat1_tmp_tmp <- t(apply(dat1_tmp, 1, function(x) (x - min(x)) / (max(x) - min(x))))
dat2_tmp_tmp <- t(apply(dat2_tmp, 1, function(x) (x - min(x)) / (max(x) - min(x))))
dat <- cbind(dat1_tmp_tmp, dat2_tmp_tmp)
dat <- na.omit(dat)

set.seed(11)
km<- kmeans(dat,2) # determin how many cluster you want, I specify 2 here
m.kmeans<- cbind(dat, km$cluster) # combine the cluster with the matrix
dim(m.kmeans)
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column

library(pheatmap)
pdf("./ZGA_zscore_kmeans/shared_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans.pdf")
pheatmap( m.kmeans[,1:ncol(m.kmeans)-1], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()

dat1_share_k4 <- dat1[rownames(m.kmeans),]
dat2_share_k27 <- dat2[rownames(m.kmeans),]

##### z-score #####
mean_value <- mean(dat1)
std_dev <- sd(as.vector(dat1))
z_score_matrix_1 <- (dat1 - mean_value) / std_dev

mean_value <- mean(dat2)
std_dev <- sd(as.vector(dat2))
z_score_matrix_2 <- (dat2 - mean_value) / std_dev

dat1_tmp2 <- scale(dat1)
dat1_tmp2[dat1_tmp2>4] <- 4
dat2_tmp2 <- scale(dat2)
dat2_tmp2[dat2_tmp2>4] <- 4

dat <- cbind(z_score_matrix_1, z_score_matrix_2)
dat[dat>4] <- 4
dat <- dat[rownames(m.kmeans),]
pdf("./ZGA_zscore_kmeans/share_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_zscore.pdf")
pheatmap( dat, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

z_score_matrix_1 <- z_score_matrix_1[rownames(m.kmeans),]
z_score_matrix_1[z_score_matrix_1 >2.5] <- 2.5
pdf("./ZGA_zscore_kmeans/share_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore.pdf")
pheatmap( z_score_matrix_1, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 

z_score_matrix_2 <- z_score_matrix_2[rownames(m.kmeans),]
z_score_matrix_2[z_score_matrix_2 >2] <- 2
pdf("./ZGA_zscore_kmeans/share_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore.pdf")
pheatmap( z_score_matrix_2, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 


############### curve ###################
dat <- data.frame(k4me3_c1=apply(z_score_matrix_1, 2, median), 
                  k27ac_c1=apply(z_score_matrix_2, 2, median), 
                  pseuo=1:69)
library(ggplot2)
pdf("./ZGA_zscore_kmeans/shared_promoter_enhancer_cocluster_curve_k4k27.pdf")
ggplot(dat, aes(x=pseuo))+
  geom_smooth(aes(y=k4me3_c1), method = "loess", color = "#9bd3b1", se = FALSE)+
  geom_smooth(aes(y=k27ac_c1), method = "loess", color = "#f4af9d", se = FALSE)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=15), #????x???̶ȱ?ǩ????????ʾ??б?Ƕ?Ϊ15?ȣ??????µ???1(hjust = 1)????????ΪTimes??СΪ20
        axis.text.y=element_text(family="Times",size=15,face="plain"), #????y???̶ȱ?ǩ???????أ???????С????????ʽΪplain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #????y????????????????
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #ȥ??Ĭ???????Ļ?ɫ??????x=0????y=0???Ӵ???ʾ(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #????ͼ?????ӱ?????????????
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #????ͼ?????ܱ?????????????
                                  size=12),
        panel.grid.major = element_blank())  #????ʾ??????
dev.off()



####################################### zscore for all ##################
all_k4 <- rbind(dat1_down_k4, dat1_up_k4, dat1_share_k4)
all_k27 <- rbind(dat2_down_k27, dat2_up_k27, dat2_share_k27)

mean_value <- mean(all_k4)
std_dev <- sd(as.vector(all_k4))
z_score_matrix_k4 <- (all_k4 - mean_value) / std_dev

mean_value <- mean(all_k27)
std_dev <- sd(as.vector(all_k27))
z_score_matrix_k27 <- (all_k27 - mean_value) / std_dev

z_score_matrix_k4[z_score_matrix_k4 >3] <- 3
pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore.pdf")
pheatmap( z_score_matrix_k4, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 


z_score_matrix_k27[z_score_matrix_k27 >3] <- 3
pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore.pdf")
pheatmap( z_score_matrix_k27, 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

# plot
pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore_1.pdf")
pheatmap( z_score_matrix_k4[1:513,], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore_2.pdf")
pheatmap( z_score_matrix_k4[514:1649,], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k4me3_zscore_3.pdf")
pheatmap( z_score_matrix_k4[1650:1808,], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#3c6aa8", "#3c6aa8","#27b4a3","#c7bd5e","#efdf31","#efdf31"))(500))
dev.off() 

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore_1.pdf")
pheatmap( z_score_matrix_k27[1:513,], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore_2.pdf")
pheatmap( z_score_matrix_k27[514:1649,], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_h3k4me3_h3k27ac_cocluster_kmeans_h3k27ac_zscore_3.pdf")
pheatmap( z_score_matrix_k27[1650:1808,], 
          cluster_rows = F, cluster_cols = F, 
          cellheight=0.2,
          show_rownames=FALSE, show_colnames=FALSE,
          color = colorRampPalette(c("#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off() 

#################################### RNA ###########################################
RNA <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/all_RNA_no_outlier.rds")
RNA_counts <- as.data.frame(RNA@assays[["RNA"]]@counts)
cds <- readRDS("/media/helab/data1/min/02_tacit/03_early_stage/20231229_allStage_integration/scRNA/RNA_trajectory.rds")
order <- sort(cds@principal_graph_aux@listData$UMAP$pseudotime)
order <- order[grep("2cell", names(order))]

RNA_counts <- RNA_counts[,grep("2cell", colnames(RNA_counts))]
RNA_counts <- RNA_counts[,names(order)]

###select genes around links
up_anno <- read.csv("./up_coclustering_ordered_anno.txt", sep = "\t")
down_anno <- read.csv("./down_coclustering_ordered_anno.txt", sep = "\t")
share_anno <- read.csv("./shared_coclustering_ordered_anno.txt", sep = "\t")
anno <- c(up_anno$Gene.Name, down_anno$Gene.Name, share_anno$Gene.Name)
RNA_counts_select <- RNA_counts[anno,]
RNA_counts_select <- apply(RNA_counts_select, 2, log10)
RNA_counts_select[RNA_counts_select == "-Inf"] <- 0

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_cocluster_RNA_3.pdf")
pheatmap(RNA_counts_select, 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = F, 
         cellheight=0.2,
         color = colorRampPalette(c("#2a628c","#3d8ec1","#f2ee75","#f9fa02"))(500))
dev.off()


pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_cocluster_RNA_1.pdf")
pheatmap(RNA_counts_select[1:374,], 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = F, 
         cellheight=0.2,
         color = colorRampPalette(c("#2a628c","#3d8ec1","#f2ee75","#f9fa02"))(500))
dev.off()

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_cocluster_RNA_2.pdf")
pheatmap(RNA_counts_select[375:1512,], 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = F, 
         cellheight=0.2,
         color = colorRampPalette(c("#2a628c","#3d8ec1","#f2ee75","#f9fa02"))(500))
dev.off()

pdf("./ZGA_zscore_kmeans/all_promoter_enhancer_cocluster_RNA_3.pdf")
pheatmap(RNA_counts_select[1513:1671,], 
         cluster_rows =F, cluster_cols = F, 
         show_rownames = F, scale ="none",
         show_colnames = F, 
         cellheight=0.2,
         color = colorRampPalette(c("#2a628c","#3d8ec1","#f2ee75","#f9fa02"))(500))
dev.off()

