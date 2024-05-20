setwd("./motifs/")

library(rvest)
library(tidyr)


##################################### calculate motif p value #####################################

html <-list.files()

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
write.table(all_motifs, "all_motifs.txt", sep = "\t", quote = F, row.names = F)