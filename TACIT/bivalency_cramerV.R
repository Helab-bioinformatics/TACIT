
########################################################################################################################
################################# 2024.4.30 bivalency score for C1 and C2 cluster ######################################
setwd("/media/helab/data3/min/TACIT/coTarget_early2cell/cicero_h3k4me3_h3k27ac/all_tss5kb/ZGA_related_cutoff0.2")

library(tidyr)
up_k4 <- read.csv("./up_promoter_enhancer_h3k4me3_2.txt", sep = "\t")
down_k4 <- read.csv("./down_promoter_enhancer_h3k4me3_2.txt", sep = "\t")
shared_k4 <- read.csv("./shared_promoter_enhancer_h3k4me3_2.txt", sep = "\t")
up_k27 <- read.csv("./up_promoter_enhancer_h3k27ac_2.txt", sep = "\t")
down_k27 <- read.csv("./down_promoter_enhancer_h3k27ac_2.txt", sep = "\t")
shared_k27 <- read.csv("./shared_promoter_enhancer_h3k27ac_2.txt", sep = "\t")

#### C1 and C2 cell id ####
pseudo_order <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/03_stage_umap/early2cell/h3k27ac_ID.txt", sep = " ")
pseudo_order <- pseudo_order[grep("22114", rownames(pseudo_order)),]
pseudo_order <- separate(pseudo_order, ID, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"), sep = "_")
pseudo_order <- unite(pseudo_order, ID2, c("V1", "V2", "V3", "V4"), sep = "_")
pseudo_order$ID <- paste0(pseudo_order$ID2, ".bam")
pseudo_order <- pseudo_order[order(pseudo_order$C1C2, pseudo_order$pseudo.time),]
C1_id <- pseudo_order$ID[which(pseudo_order$C1C2 == "C1")]
C2_id <- pseudo_order$ID[which(pseudo_order$C1C2 == "C2")]

############# h3k4me3 and h3k27ac raw signals in C1 and C2 regions ####################
################# down ########################
down_k4 <- down_k4[,2:90]
down_k27 <- down_k27[,2:90]
down_k4_c1 <- down_k4[1:58,]; down_k4_c2 <- down_k4[59:89,]
down_k27_c1 <- down_k27[1:58,]; down_k27_c2 <- down_k27[59:89,]

down_k4_c1_signal <- apply(down_k4_c1, 1, mean)
down_k4_c2_signal <- apply(down_k4_c2, 1, mean)
down_k27_c1_signal <- apply(down_k27_c1, 1, mean)
down_k27_c2_signal <- apply(down_k27_c2, 1, mean)
dat <- data.frame(k4_signal=c(down_k4_c1_signal, down_k4_c2_signal), 
                  k27_signal=c(down_k27_c1_signal, down_k27_c2_signal), 
                  celltype=c(rep("c1", 58), rep("c2", 31)))
library(tidyr)
dat <- gather(dat, marker, value, -celltype)

library(ggplot2)
library(ggunchained)
pdf("./ZGA_kmeans/down_peaks_k4_k27ac_counts_average.pdf", width = 10)
ggplot(dat,aes(x=celltype,y=value,fill=marker))+geom_split_violin()+
  ylim(0,1)+
  scale_fill_manual(values = c("#7697CB","#D293C1"))+
  theme_bw()
dev.off()

################# shared ########################
shared_k4 <- shared_k4[,2:90]
shared_k27 <- shared_k27[,2:90]
shared_k4_c1 <- shared_k4[1:58,]; shared_k4_c2 <- shared_k4[59:89,]
shared_k27_c1 <- shared_k27[1:58,]; shared_k27_c2 <- shared_k27[59:89,]

shared_k4_c1_signal <- apply(shared_k4_c1, 1, mean)
shared_k4_c2_signal <- apply(shared_k4_c2, 1, mean)
shared_k27_c1_signal <- apply(shared_k27_c1, 1, mean)
shared_k27_c2_signal <- apply(shared_k27_c2, 1, mean)
dat <- data.frame(k4_signal=c(shared_k4_c1_signal, shared_k4_c2_signal), 
                  k27_signal=c(shared_k27_c1_signal, shared_k27_c2_signal), 
                  celltype=c(rep("c1", 58), rep("c2", 31)))
library(tidyr)
dat <- gather(dat, marker, value, -celltype)

library(ggplot2)
#library(ggpol)
pdf("./ZGA_kmeans/shared_peaks_k4_k27ac_counts_average.pdf", width = 10)
ggplot(dat,aes(x=celltype,y=value,fill=marker))+geom_split_violin()+
  scale_fill_manual(values = c("#7697CB","#D293C1"))+
  ylim(0,1.2)+
  theme_bw()
dev.off()


################# up ########################
up_k4 <- up_k4[,2:90]
up_k27 <- up_k27[,2:90]
up_k4_c1 <- up_k4[1:58,]; up_k4_c2 <- up_k4[59:89,]
up_k27_c1 <- up_k27[1:58,]; up_k27_c2 <- up_k27[59:89,]

up_k4_c1_signal <- apply(up_k4_c1, 1, mean)
up_k4_c2_signal <- apply(up_k4_c2, 1, mean)
up_k27_c1_signal <- apply(up_k27_c1, 1, mean)
up_k27_c2_signal <- apply(up_k27_c2, 1, mean)
dat <- data.frame(k4_signal=c(up_k4_c1_signal, up_k4_c2_signal), 
                  k27_signal=c(up_k27_c1_signal, up_k27_c2_signal), 
                  celltype=c(rep("c1", 58), rep("c2", 31)))
library(tidyr)
dat <- gather(dat, marker, value, -celltype)

library(ggplot2)
pdf("./ZGA_kmeans/up_peaks_k4_k27ac_counts_average.pdf", width = 10)
ggplot(dat,aes(x=celltype,y=value,fill=marker))+geom_split_violin()+
  ylim(0,0.9)+
  scale_fill_manual(values = c("#7697CB","#D293C1"))+
  theme_bw()
dev.off()


########################################## bivalency score ###########################################

##### down ####
dim <- nrow(down_k4) * ncol(down_k4)
down_biva <- matrix(1:dim, ncol=ncol(down_k4))
down_biva <- as.data.frame(down_biva)
for (i in 1:ncol(down_k4)) {
  sum <- down_k4[,i] + down_k27[,i]
  dif <- abs(down_k4[,i] - down_k27[,i]) + 1
  down_biva[,i] <- (sum / dif)
}

## 区分C1和C2
down_biva <- apply(down_biva, 2, as.numeric)
down_biva_C1 <- c(down_biva[,1:58])# C1细胞是58个
down_biva_C2 <- c(down_biva[,59:89])
down_biva_C1 <- down_biva_C1[down_biva_C1<2]
down_biva_C2 <- down_biva_C2[down_biva_C2<2]
down_biva_all <- data.frame(bivalency=c(down_biva_C1, down_biva_C2),
                            groshared=c(rep("c1", length(down_biva_C1)), rep("c2", length(down_biva_C2))))
library(ggplot2)
library(ggpubr)
pdf("./ZGA_kmeans/bivalency_down.pdf")
ggplot(data=down_biva_all,mapping = aes(x = groshared, y = bivalency))+
  geom_violin(aes(fill = groshared), trim = T) + 
  geom_boxplot(width = 0.05)+
  theme(legend.position = "none")+
  stat_compare_means()+
  theme_bw()
dev.off()

## heatmap ##
id <- seq(from=1,to=85,by=5) # length(id) = 17
dim <- nrow(down_biva) * 17
down_biva_smo <- matrix(1:dim,ncol=17)
n=1
for (i in id) {
  m=i+4
  down_biva_smo[,n] <- rowSums(down_biva[,i:m])
  n=n+1
}
pdf("./ZGA_kmeans/bivalency_down_smoth.pdf")
pheatmap(down_biva_smo, 
         cluster_rows =T, cluster_cols = F, 
         show_rownames = T, scale ="row",
         show_colnames = F, 
         cellheight=0.2,
         color = colorRampPalette(c("#0561c6","#0561c6","#f4e4b8","#e2904b","#ac2728"))(500))
dev.off()


################# shared ########################
shared_k4 <- shared_k4[,2:90]
shared_k27 <- shared_k27[,2:90]

## bivalency score
dim <- nrow(shared_k4) * ncol(shared_k4)
shared_biva <- matrix(1:dim, ncol=ncol(shared_k4))
shared_biva <- as.data.frame(shared_biva)
for (i in 1:ncol(shared_k4)) {
  sum <- shared_k4[,i] + shared_k27[,i]
  dif <- abs(shared_k4[,i] - shared_k27[,i]) + 1
  shared_biva[,i] <- (sum / dif)
}

## C1
shared_biva <- apply(shared_biva, 2, as.numeric)
shared_biva_C1 <- c(shared_biva[,1:58])# C1细胞是58个
shared_biva_C2 <- c(shared_biva[,59:89])
shared_biva_C1 <- shared_biva_C1[shared_biva_C1<2]
shared_biva_C2 <- shared_biva_C2[shared_biva_C2<2]
shared_biva_all <- data.frame(bivalency=c(shared_biva_C1, shared_biva_C2),
                          groshared=c(rep("c1", length(shared_biva_C1)), rep("c2", length(shared_biva_C2))))
library(ggplot2)
library(ggpubr)
pdf("./ZGA_kmeans/bivalency_shared.pdf")
ggplot(data=shared_biva_all,mapping = aes(x = groshared, y = bivalency))+
  geom_violin(aes(fill = groshared), trim = T) + 
  geom_boxplot(width = 0.05)+
  theme(legend.position = "none")+
  stat_compare_means()+
  theme_bw()
dev.off()


##### 对于每个peak取在C1细胞群或者C2细胞群中的median的bivalency score
##### down ####
dim <- nrow(down_k4) * ncol(down_k4)
down_biva <- matrix(1:dim, ncol=ncol(down_k4))
down_biva <- as.data.frame(down_biva)
for (i in 1:ncol(down_k4)) {
  sum <- down_k4[,i] + down_k27[,i]
  dif <- abs(down_k4[,i] - down_k27[,i]) + 1
  down_biva[,i] <- (sum / dif)
}

down_biva_C1 <- down_biva[,1:58]# C1细胞是58个
down_biva_C2 <- down_biva[,59:89]
down_biva_C1 <- apply(down_biva_C1, 1, median)
down_biva_C2 <- apply(down_biva_C2, 1, median)
dat <- data.frame(bivalency=c(down_biva_C1, down_biva_C2), celltype=rep(c("C1", "C2"), each=515))

library(ggplot2)
library(ggpubr)
pdf("./ZGA_kmeans/bivalency_down_average.pdf")
ggplot(data=dat,mapping = aes(x = celltype, y = bivalency))+
  geom_violin(aes(fill = celltype), trim = T) + 
  geom_boxplot(width = 0.05)+
  theme(legend.position = "none")+
  stat_compare_means()+
  theme_bw()
dev.off()

#######################################################################################
##################################### binarization and cramerV similarity #############

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
    Chi.sq = ssharedpressWarnings(chisq.test(x, y, correct=FALSE, ...)$statistic)
    Phi    = Chi.sq / N
    Row    = length(unique(x))
    C      = length(unique(y))
    CV     =  sqrt(Phi / min(Row-1, C-1))
  }
  
  if(is.matrix(x)){x=as.table(x)}
  if(is.table(x)){
    TABLE  = x
    N      = sum(TABLE)
    Chi.sq = ssharedpressWarnings(chisq.test(TABLE, correct=FALSE, ...)$statistic)
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
    return(data.frame(Cramer.V=CV, lower.ci=NA, sharedper.ci=NA))} 
  
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
        Chi.sq = ssharedpressWarnings(chisq.test(Input[,1], Input[,2], 
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
  if(ci==TRUE){return(data.frame(Cramer.V=CV, lower.ci=CI1, sharedper.ci=CI2))}  
}

##### down ####
down_k4[down_k4 > 0] <- 1
down_k4[down_k4 == 0] <- 0
down_k27[down_k27 > 0] <- 1
down_k27[down_k27 == 0] <- 2

cramer_c1 <- c()
for (i in 1:58) {
  cramer_c1[i] <- cramerV(down_k4[,i], down_k27[,i])
}

cramer_c2 <- c()
for (i in 59:89) {
  cramer_c2[i] <- cramerV(down_k4[,i], down_k27[,i])
}
cramer_c2 <- na.omit(cramer_c2)

dat <- data.frame(cramerv=c(cramer_c1, cramer_c2), celltype=c(rep("C1",58), rep("C2",31)))
pdf("./ZGA_kmeans/cramerV_down.pdf")
ggplot(data=dat,mapping = aes(x = celltype, y = cramerv))+
  geom_violin(aes(fill = celltype), trim = T) + 
  geom_boxplot(width = 0.05)+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test")+
  theme_bw()
dev.off()


##### shared ####
shared_k4[shared_k4 > 0] <- 1
shared_k4[shared_k4 == 0] <- 0
shared_k27[shared_k27 > 0] <- 1
shared_k27[shared_k27 == 0] <- 2

cramer_c1 <- c()
for (i in 1:58) {
  cramer_c1[i] <- cramerV(shared_k4[,i], shared_k27[,i])
}

cramer_c2 <- c()
for (i in 59:89) {
  cramer_c2[i] <- cramerV(shared_k4[,i], shared_k27[,i])
}
cramer_c2 <- na.omit(cramer_c2)

dat <- data.frame(cramerv=c(cramer_c1, cramer_c2), celltype=c(rep("C1",58), rep("C2",31)))
pdf("./ZGA_kmeans/cramerV_shared.pdf")
ggplot(data=dat,mapping = aes(x = celltype, y = cramerv))+
  geom_violin(aes(fill = celltype), trim = T) + 
  geom_boxplot(width = 0.05)+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test")+
  theme_bw()
dev.off()



##### shared ####
shared_k4 <- shared_k4[,2:90]
shared_k27 <- shared_k27[,2:90]

shared_k4[shared_k4 > 0] <- 1
shared_k4[shared_k4 == 0] <- 0
shared_k27[shared_k27 > 0] <- 1
shared_k27[shared_k27 == 0] <- 2

cramer_c1 <- c()
for (i in 1:58) {
  cramer_c1[i] <- cramerV(shared_k4[,i], shared_k27[,i])
}

cramer_c2 <- c()
for (i in 59:89) {
  cramer_c2[i] <- cramerV(shared_k4[,i], shared_k27[,i])
}
cramer_c2 <- na.omit(cramer_c2)

dat <- data.frame(cramerv=c(cramer_c1, cramer_c2), celltype=c(rep("C1",58), rep("C2",31)))
pdf("./ZGA_kmeans/cramerV_shared.pdf")
ggplot(data=dat,mapping = aes(x = celltype, y = cramerv))+
  geom_violin(aes(fill = celltype), trim = T) + 
  geom_boxplot(width = 0.05)+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test")+
  theme_bw()
dev.off()


###################################### binarization and calculate percentage ##############################
up_k4 <- read.csv("./up_promoter_enhancer_h3k4me3.txt", sep = "\t")
down_k4 <- read.csv("./down_promoter_enhancer_h3k4me3.txt", sep = "\t")
shared_k4 <- read.csv("./shared_promoter_enhancer_h3k4me3.txt", sep = "\t")
up_k27 <- read.csv("./up_promoter_enhancer_h3k27ac.txt", sep = "\t")
down_k27 <- read.csv("./down_promoter_enhancer_h3k27ac.txt", sep = "\t")
shared_k27 <- read.csv("./shared_promoter_enhancer_h3k27ac.txt", sep = "\t")

#### C1 and C2 cell id ####
pseudo_order <- read.csv("/media/helab/data1/min/02_tacit/03_early_stage/all_h3k27ac/03_stage_umap/early2cell/h3k27ac_ID.txt", sep = " ")
pseudo_order <- pseudo_order[grep("22114", rownames(pseudo_order)),]
pseudo_order <- separate(pseudo_order, ID, c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10"), sep = "_")
pseudo_order <- unite(pseudo_order, ID2, c("V1", "V2", "V3", "V4"), sep = "_")
pseudo_order$ID <- paste0(pseudo_order$ID2, ".bam")
pseudo_order <- pseudo_order[order(pseudo_order$C1C2, pseudo_order$pseudo.time),]
C1_id <- pseudo_order$ID[which(pseudo_order$C1C2 == "C1")]
C2_id <- pseudo_order$ID[which(pseudo_order$C1C2 == "C2")]

################################# down  (C1) ########################################
down_k27_C1 <- as.matrix(down_k27[,which(colnames(down_k27) %in% C1_id)])
down_k4_C1 <- as.matrix(down_k4[,which(colnames(down_k4) %in% C1_id)])
# 计算每一行中两个矩阵在相同位置都为1的数量
both_one <- rowSums(down_k4_C1 == 1 & down_k27_C1 == 1)
# 计算每一行中两个矩阵在相同位置一个为1一个为0的数量
one_zero <- rowSums((down_k4_C1 == 1 & down_k27_C1 == 0) | (down_k4_C1 == 0 & down_k27_C1 == 1))
# 计算每一行中两个矩阵在相同位置都为0的数量
both_zero <- rowSums(down_k4_C1 == 0 & down_k27_C1 == 0)
# 计算总比例
total_positions <- ncol(down_k4_C1)  # 总位置数，这里是每行的元素数量，即89
# 计算比例
proportion_both_one <- both_one / total_positions
proportion_one_zero <- one_zero / total_positions
proportion_both_zero <- both_zero / total_positions

# 查看结果
results_down_c1 <- data.frame(BothOne = proportion_both_one, OneZero = proportion_one_zero, BothZero = proportion_both_zero)

### down  (c2) ####
down_k27_c2 <- as.matrix(down_k27[,which(colnames(down_k27) %in% C2_id)])
down_k4_c2 <- as.matrix(down_k4[,which(colnames(down_k4) %in% C2_id)])
# 计算每一行中两个矩阵在相同位置都为1的数量
both_one <- rowSums(down_k4_c2 == 1 & down_k27_c2 == 1)
# 计算每一行中两个矩阵在相同位置一个为1一个为0的数量
one_zero <- rowSums((down_k4_c2 == 1 & down_k27_c2 == 0) | (down_k4_c2 == 0 & down_k27_c2 == 1))
# 计算每一行中两个矩阵在相同位置都为0的数量
both_zero <- rowSums(down_k4_c2 == 0 & down_k27_c2 == 0)
# 计算总比例
total_positions <- ncol(down_k4_c2)  # 总位置数，这里是每行的元素数量，即89
# 计算比例
proportion_both_one <- both_one / total_positions
proportion_one_zero <- one_zero / total_positions
proportion_both_zero <- both_zero / total_positions

# 查看结果
results_down_c2 <- data.frame(BothOne = proportion_both_one, OneZero = proportion_one_zero, BothZero = proportion_both_zero)

bothone_down <- data.frame(bothone_c1=results_down_c1$BothOne, bothone_c2=results_down_c2$BothOne)

###################################### up  (C1) #############################################
up_k27_C1 <- as.matrix(up_k27[,which(colnames(up_k27) %in% C1_id)])
up_k4_C1 <- as.matrix(up_k4[,which(colnames(up_k4) %in% C1_id)])
# 计算每一行中两个矩阵在相同位置都为1的数量
both_one <- rowSums(up_k4_C1 == 1 & up_k27_C1 == 1)
# 计算每一行中两个矩阵在相同位置一个为1一个为0的数量
one_zero <- rowSums((up_k4_C1 == 1 & up_k27_C1 == 0) | (up_k4_C1 == 0 & up_k27_C1 == 1))
# 计算每一行中两个矩阵在相同位置都为0的数量
both_zero <- rowSums(up_k4_C1 == 0 & up_k27_C1 == 0)
# 计算总比例
total_positions <- ncol(up_k4_C1)  # 总位置数，这里是每行的元素数量，即89
# 计算比例
proportion_both_one <- both_one / total_positions
proportion_one_zero <- one_zero / total_positions
proportion_both_zero <- both_zero / total_positions

# 查看结果
results_up_c1 <- data.frame(BothOne = proportion_both_one, OneZero = proportion_one_zero, BothZero = proportion_both_zero)

### up  (c2) ####
up_k27_c2 <- as.matrix(up_k27[,which(colnames(up_k27) %in% C2_id)])
up_k4_c2 <- as.matrix(up_k4[,which(colnames(up_k4) %in% C2_id)])
# 计算每一行中两个矩阵在相同位置都为1的数量
both_one <- rowSums(up_k4_c2 == 1 & up_k27_c2 == 1)
# 计算每一行中两个矩阵在相同位置一个为1一个为0的数量
one_zero <- rowSums((up_k4_c2 == 1 & up_k27_c2 == 0) | (up_k4_c2 == 0 & up_k27_c2 == 1))
# 计算每一行中两个矩阵在相同位置都为0的数量
both_zero <- rowSums(up_k4_c2 == 0 & up_k27_c2 == 0)
# 计算总比例
total_positions <- ncol(up_k4_c2)  # 总位置数，这里是每行的元素数量，即89
# 计算比例
proportion_both_one <- both_one / total_positions
proportion_one_zero <- one_zero / total_positions
proportion_both_zero <- both_zero / total_positions

# 查看结果
results_up_c2 <- data.frame(BothOne = proportion_both_one, OneZero = proportion_one_zero, BothZero = proportion_both_zero)

bothone_up <- data.frame(bothone_c1=results_up_c1$BothOne, bothone_c2=results_up_c2$BothOne)



###################################### shared  (C1) #############################################
shared_k27_C1 <- as.matrix(shared_k27[,which(colnames(shared_k27) %in% C1_id)])
shared_k4_C1 <- as.matrix(shared_k4[,which(colnames(shared_k4) %in% C1_id)])
# 计算每一行中两个矩阵在相同位置都为1的数量
both_one <- rowSums(shared_k4_C1 == 1 & shared_k27_C1 == 1)
# 计算每一行中两个矩阵在相同位置一个为1一个为0的数量
one_zero <- rowSums((shared_k4_C1 == 1 & shared_k27_C1 == 0) | (shared_k4_C1 == 0 & shared_k27_C1 == 1))
# 计算每一行中两个矩阵在相同位置都为0的数量
both_zero <- rowSums(shared_k4_C1 == 0 & shared_k27_C1 == 0)
# 计算总比例
total_positions <- ncol(shared_k4_C1)  # 总位置数，这里是每行的元素数量，即89
# 计算比例
proportion_both_one <- both_one / total_positions
proportion_one_zero <- one_zero / total_positions
proportion_both_zero <- both_zero / total_positions

# 查看结果
results_shared_c1 <- data.frame(BothOne = proportion_both_one, OneZero = proportion_one_zero, BothZero = proportion_both_zero)

### shared  (c2) ####
shared_k27_c2 <- as.matrix(shared_k27[,which(colnames(shared_k27) %in% C2_id)])
shared_k4_c2 <- as.matrix(shared_k4[,which(colnames(shared_k4) %in% C2_id)])
# 计算每一行中两个矩阵在相同位置都为1的数量
both_one <- rowSums(shared_k4_c2 == 1 & shared_k27_c2 == 1)
# 计算每一行中两个矩阵在相同位置一个为1一个为0的数量
one_zero <- rowSums((shared_k4_c2 == 1 & shared_k27_c2 == 0) | (shared_k4_c2 == 0 & shared_k27_c2 == 1))
# 计算每一行中两个矩阵在相同位置都为0的数量
both_zero <- rowSums(shared_k4_c2 == 0 & shared_k27_c2 == 0)
# 计算总比例
total_positions <- ncol(shared_k4_c2)  # 总位置数，这里是每行的元素数量，即89
# 计算比例
proportion_both_one <- both_one / total_positions
proportion_one_zero <- one_zero / total_positions
proportion_both_zero <- both_zero / total_positions

# 查看结果
results_shared_c2 <- data.frame(BothOne = proportion_both_one, OneZero = proportion_one_zero, BothZero = proportion_both_zero)

bothone_shared <- data.frame(bothone_c1=results_shared_c1$BothOne, bothone_c2=results_shared_c2$BothOne)


########################################## plot #####################################
bothone_up$group <- "up"
bothone_down$group <- "down"
bothone_shared$group <- "shared"
dat <- rbind(bothone_up, bothone_down, bothone_shared)

library(tidyr)
dat <- gather(dat, celltype, value, -group)

library(ggplot2)
pdf("./ZGA_kmeans/up_down_shared_bothK27K4_portion.pdf", width = 10)
ggplot(dat,aes(x=group,y=value,fill=celltype))+geom_split_violin()+
  scale_fill_manual(values = c("#7697CB","#D293C1"))+
  theme_bw()
dev.off()

############ 连线图 #########
library(ggplot2)
# 假设你的数据框是df，有ID列标识每行
df <- data.frame(ID = paste0("loci", 1:nrow(bothone_up)), A = bothone_up$bothone_c1, B = bothone_up$bothone_c2)
# 将数据从宽格式转换为长格式
df_long <- reshape2::melt(df, id.vars = "ID", measure.vars = c("A", "B"))
pdf("./ZGA_kmeans/up_bothk27k4_portion.pdf")
ggplot(df_long, aes(x = variable, y = value, group = ID)) +
  geom_line(color = "lightgray") +  # 连接线
  geom_point() +  # 可选，增加点
  theme_bw() +  # 选择一个简洁的主题
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1.2) +  # 平均值连线
  stat_compare_means(comparisons = list(c("A", "B"))) +  # 统计显著性
  labs(x = "组别", y = "数值大小", title = "AB两列数值的连接线图")
dev.off()



