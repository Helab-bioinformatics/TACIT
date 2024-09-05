setwd("./04_ICM_TE_active/TE_only_active")
webshot::install_phantomjs()
library(networkD3)
library(tidyr)

###binned bed file should be sorted at first
A <- read.csv("./4cell_10_chrhmm_inte_TE_only_active.txt", header = F, sep = "\t")
B <- read.csv("./8cell_10_chrhmm_inte_TE_only_active.txt", header = F, sep = "\t")
C <- read.csv("./morula_10_chrhmm_inte_TE_only_active.txt", header = F, sep = "\t")
D <- read.csv("./TE_10_chrhmm_inte_TE_only_active.txt", header = F, sep = "\t")

A[which(A$V4 == "E1"),5] <- "Pa"
A[which(A$V4 == "E2"),5] <- "Es"
A[which(A$V4 == "E3"),5] <- "Un"
A[which(A$V4 == "E4"),5] <- "Es"
A[which(A$V4 == "E5"),5] <- "Ep"
A[which(A$V4 == "E6"),5] <- "Ti"
A[which(A$V4 == "E7"),5] <- "Te"
A[which(A$V4 == "E8"),5] <- "Un"
A[which(A$V4 == "E9"),5] <- "K27"
A[which(A$V4 == "E10"),5] <- "Un"

B[which(B$V4 == "E1"),5] <- "Pa"
B[which(B$V4 == "E2"),5] <- "Es"
B[which(B$V4 == "E3"),5] <- "Un"
B[which(B$V4 == "E4"),5] <- "Es"
B[which(B$V4 == "E5"),5] <- "Ep"
B[which(B$V4 == "E6"),5] <- "Ti"
B[which(B$V4 == "E7"),5] <- "Te"
B[which(B$V4 == "E8"),5] <- "Un"
B[which(B$V4 == "E9"),5] <- "K27"
B[which(B$V4 == "E10"),5] <- "Un"

C[which(C$V4 == "E1"),5] <- "Pa"
C[which(C$V4 == "E2"),5] <- "Es"
C[which(C$V4 == "E3"),5] <- "Un"
C[which(C$V4 == "E4"),5] <- "Es"
C[which(C$V4 == "E5"),5] <- "Ep"
C[which(C$V4 == "E6"),5] <- "Ti"
C[which(C$V4 == "E7"),5] <- "Te"
C[which(C$V4 == "E8"),5] <- "Un"
C[which(C$V4 == "E9"),5] <- "K27"
C[which(C$V4 == "E10"),5] <- "Un"

D[which(D$V4 == "E1"),5] <- "Pa"
D[which(D$V4 == "E2"),5] <- "Es"
D[which(D$V4 == "E3"),5] <- "Un"
D[which(D$V4 == "E4"),5] <- "Es"
D[which(D$V4 == "E5"),5] <- "Ep"
D[which(D$V4 == "E6"),5] <- "Ti"
D[which(D$V4 == "E7"),5] <- "Te"
D[which(D$V4 == "E8"),5] <- "Un"
D[which(D$V4 == "E9"),5] <- "K27"
D[which(D$V4 == "E10"),5] <- "Un"

A <- unite(A, "min", c("V1", "V2"), sep = ":", remove = T) 
A <- unite(A, "loci", c("min", "V3"), sep = "-", remove = T) 
B <- unite(B, "min", c("V1", "V2"), sep = ":", remove = T) 
B <- unite(B, "loci", c("min", "V3"), sep = "-", remove = T) 
C <- unite(C, "min", c("V1", "V2"), sep = ":", remove = T) 
C <- unite(C, "loci", c("min", "V3"), sep = "-", remove = T) 
D <- unite(D, "min", c("V1", "V2"), sep = ":", remove = T) 
D <- unite(D, "loci", c("min", "V3"), sep = "-", remove = T) 

identical(A$loci, B$loci)
identical(B$loci, C$loci)
identical(C$loci, D$loci)
# four cell type ===========================================================================
dat <- cbind(A,B,C,D)
dat <- dat[,c(1, 3, 6, 9, 12)]
# colnames(dat) <- c("loci","4cell", "8cell", "morula", "TE")
# write.csv(dat, "link_direction_4cell_TE.csv", quote = F)

colnames(dat) <- c("loci","L1", "L2", "L3", "L4")
dat$L1 <- paste0("1", dat$L1)
dat$L2 <- paste0("2", dat$L2)
dat$L3 <- paste0("3", dat$L3)
dat$L4 <- paste0("4", dat$L4)
dat$value <- 1
level1 <- aggregate(dat$value, by=list(dat$L1,dat$L2), sum)
level2 <- aggregate(dat$value, by=list(dat$L2,dat$L3), sum)
level3 <- aggregate(dat$value, by=list(dat$L3,dat$L4), sum)

links <- rbind(level3, level2, level1)
colnames(links) <- c('source',"target","value")
tail(links)

nodes <- data.frame(name=unique(c(links$source,links$target)))
links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1

Sankey.p <- sankeyNetwork(Links = links, Nodes = nodes,
                          Source = "IDsource", Target = "IDtarget",
                          Value = "value", NodeID = "name", fontFamily = "Arial",
                          fontSize = 12, nodeWidth = 50, nodePadding = 5,
                          height = 200, width = 500,
                          sinksRight=F)
Sankey.p

saveNetwork(Sankey.p, "./4cellToTE_TE_only_active_3.html", selfcontained = TRUE)
webshot("4cellToTE_TE_only_active_3.html" , "./4cellToTE_TE_only_active_3.pdf")


install.packages("webshot",dependencies = TRUE)
library(webshot)
webshot::install_phantomjs()
