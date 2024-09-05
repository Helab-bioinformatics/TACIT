
##################################### processing results from Homer #####################################
library(rvest)
library(tidyr)

html <- list.files()
for (i in 1:length(html)) {
  page <- read_html(paste0(html[i], "/homerResults.html"),)
  data <- page %>%
    html_node("table") %>%
    html_table()
  data <- separate(data, X8, c("TF", "xx"), sep = "/")
  data <- separate(data, X1, c("X1", "na"), sep = "\n")
  write.table(data[,c("X1","X3", "TF")], paste0(html[i], "_homer.txt"), sep = "\t", quote = F, row.names = F)
  
  links <- page %>% html_nodes("a") %>% html_attr("href")
  links <- links[grep(".info.html", links)]
  similarTF <- data.frame(TF=c("min"), level=c("min"))
  for (j in 1:length(links)) {
    path <- getwd()
    page2 <- read_html(paste0(path, "/",html[i], "/", links[j]))
    data2 <- page2 %>%
      html_table()
    data2 <- as.data.frame(data2[2])
    num <- nrow(data2)/7
    data2 <- data2[seq(1, nrow(data2), by=7),]
    data2 <- separate(data2, X1, c("TF", "xx"), sep = "/")
    mat <- data.frame(TF=data2$TF, level=paste0("motif", j))
    similarTF <- rbind(similarTF, mat)
  }
  write.table(similarTF[-1,], paste0(html[i], "_homer_similar_TFs.txt"), sep = "\t", quote = F, row.names = F)
}

mat <- list()
for (i in 1:length(html)) {
  a <- read.csv(paste0(html[i], "_homer.txt"), sep = "\t", header = T)
  colnames(a) <- a[1,]
  a <- a[-1,]
  a$Rank <- paste0("motif", a$Rank)
  a_similar <- read.csv(paste0(html[i], "_homer_similar_TFs.txt"), sep = "\t", header = T)
  a_similar$pvalue <- a[match(a_similar$level, a$Rank),2]
  a_similar$group <- paste0(html[i])
  a_similar <- a_similar[!duplicated(a_similar$TF),]
  mat[[i]] <- a_similar
}
all_motifs <- do.call(rbind, mat)
write.table(all_motifs, "c8_1_motifs.txt", sep = "\t", quote = F, row.names = F)

########################################## plot ########################################
library(tidyr)
tf_plot <- read.csv("./tf_plot.txt", sep = "\t")
rownames(tf_plot) <- tf_plot[,1]
tf_plot <- tf_plot[,-1]
plot_tf <- tf_plot[,1:7]
plot_tf$gene <- rownames(plot_tf)
plot_tf <- gather(plot_tf, key=cell, value=pvalue, -c(gene))
plot_tf <- unite(plot_tf, id, c("cell", "gene"), sep = "_", remove = F)

plot_exp <- tf_plot[,9:15]
colnames(plot_exp) <- colnames(tf_plot)[1:7]
plot_exp$gene <- rownames(plot_exp)
plot_exp <- gather(plot_exp, key=cell, value=exp, -c(gene))
plot_exp <- unite(plot_exp, id, c("cell", "gene"), sep = "_", remove = F)
plot_exp$pvalue <- plot_tf$pvalue[match(plot_exp$id, plot_tf$id)]

plot_exp$gene <- factor(plot_exp$gene, levels = c(rownames(tf_plot)))
plot_exp$cell <- factor(plot_exp$cell, levels = rev(c("C2_1", "C2_2","C4_2",  "C4_4", "C4_1", "C4_3", "C8_all")))
plot_exp[which(plot_exp$pvalue >15),5] <- 15

plot_exp[which(plot_exp$exp == 0),6] <- "L1"
plot_exp[which(plot_exp$exp >0 & plot_exp$exp<10),6] <- "L2"
plot_exp[which(plot_exp$exp >10 & plot_exp$exp<50),6] <- "L3"
plot_exp[which(plot_exp$exp >50 & plot_exp$exp<100),6] <- "L4"
plot_exp[which(plot_exp$exp >100 ),6] <- "L5"
colnames(plot_exp) <- c(colnames(plot_exp[,1:5]), "group")

plot_exp$logexp <- log10(plot_exp$exp + 1)
plot_exp[which(plot_exp$logexp > 2.0),7] <- 2.0

library(ggplot2)
pdf("./previous_top1957_motifs_dot_2.pdf", width = 12, height = 4)
ggplot(plot_exp,aes(x=gene,y=cell))+
  geom_point(aes(size=`pvalue`,
                 color=`group`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  labs(x=NULL,y=NULL)
dev.off()

