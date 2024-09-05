setwd("./cotacit_early2cell/cramerV/")

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

####################################### for sub regions ###################################
bed <- list.files("./cotacit_early2cell/cramerV")
bed <- bed[grep("2kb.mm10.bed", bed)]
bed <- bed[grep("pdf",invert = T, bed)]

for ( j in 1:7) {
  library(cisTopic)
  pathToBams1 <- c('../01_h3k4me3/rename/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("22114",bamFiles1)]
  regions <- bed[j]
  cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac', min.cells = 0, min.regions = 0)
  h3k4me3 <- as.matrix(cisTopicObject@count.matrix)
  h3k4me3[h3k4me3 > 1] <- 2
  
  pathToBams1 <- c('../02_h3k27ac/rename/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("22114",bamFiles1)]
  cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac', min.cells = 0, min.regions = 0)
  h3k27ac <- as.matrix(cisTopicObject@count.matrix)
  h3k27ac[h3k27ac > 1] <- 2
  
  pathToBams1 <- c('../03_h3k27me3/rename/')
  bamFiles1 <- paste(pathToBams1, list.files(pathToBams1), sep='')
  bamFiles1 <- bamFiles1[grep("22114",bamFiles1)]
  cisTopicObject <- createcisTopicObjectFromBAM(bamFiles1, regions, paired = T, project.name='E2.5_h3k27ac', min.cells = 0, min.regions = 0)
  h3k27me3 <- as.matrix(cisTopicObject@count.matrix)
  h3k27me3[h3k27me3 > 1] <- 2
  
  cramer_k4k27ac <- c()
  for (i in 1:89) {
    cramer_k4k27ac[i] <- cramerV(h3k4me3[,i], h3k27ac[,i])
  }
  names(cramer_k4k27ac) <- colnames(h3k4me3)
  
  cramer_k4k27me3 <- c()
  for (i in 1:89) {
    cramer_k4k27me3[i] <- cramerV(h3k4me3[,i], h3k27me3[,i])
  }
  names(cramer_k4k27me3) <- colnames(h3k4me3)
  
  cramer_k27k27 <- c()
  for (i in 1:89) {
    cramer_k27k27[i] <- cramerV(h3k27ac[,i], h3k27me3[,i])
  }
  names(cramer_k27k27) <- colnames(h3k27ac)
  
  dat <- data.frame(k4k27ac=cramer_k4k27ac, k4k27me3=cramer_k4k27me3, k27ack27me3=cramer_k27k27)
  dat$id <- rownames(dat)
  library(tidyr)
  dat <- gather(dat, category, value, k4k27ac, k4k27me3, k27ack27me3)
  dat <- dat[which(dat$value <0.03),] #filtering
  
  library(ggplot2)
  library(ggpubr)
  pdf(paste0(bed[j], "_cramerV_2kbbin.pdf"))
  ggplot(data=dat,mapping = aes(x = category, y = value))+
    geom_violin(aes(fill = category), trim = T) + 
    geom_boxplot(width = 0.05)+
    ylim(0, 0.03)+
    stat_compare_means()+
    scale_fill_manual(values = c("#e8490f","#f18800","#e4ce00","#9ec417","#13a983","#44c1f0"))+
    theme(legend.position = "none")+
    theme_bw()
  dev.off()
  
}



