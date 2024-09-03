
#### select region for calculating cramerV similarity###
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


### rough classification
all_anno <- read.csv("./top780_ICM_TE_anno.txt", sep = "\t")

dim <- ncol(all_anno) * ncol(all_anno)
cramer_matrix <- matrix(1:dim,ncol = ncol(all_anno))
for (i in 1:ncol(all_anno)) {
  for (j in i:ncol(all_anno)) {
    cramer_matrix[i,j] <- cramerV(all_anno[,i], all_anno[,j])
    cramer_matrix[j,i] <- cramer_matrix[i,j]
  }
}
rownames(cramer_matrix) <- colnames(all_anno)
colnames(cramer_matrix) <- colnames(all_anno)

cramer_matrix_ICM <- cramer_matrix[grep("ICM", rownames(cramer_matrix)),]
cramer_matrix_TE <- cramer_matrix[grep("TE", rownames(cramer_matrix)),]

carmer_ICM_ave <- apply(cramer_matrix_ICM, 2, mean)
carmer_TE_ave <- apply(cramer_matrix_TE, 2, mean)

library(ggplot2)
dat <- data.frame(cells=names(carmer_ICM_ave),icm_score=carmer_ICM_ave, te_score=carmer_TE_ave)

pdf("allStage_cramerV_ICM_TE.pdf")
ggplot(data=dat, aes(x=icm_score,y=te_score)) +
  geom_point()  + 
  xlim(0.01,0.6)+
  ylim(0.01,0.6)+
  geom_text(aes(label = cells)) +
  theme_bw()
dev.off()
