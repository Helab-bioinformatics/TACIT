
############################## enhancer strong ################################
rm(list = ls())
stage <- c("zygote",  paste0("2cell_", 1:2), paste0("4cell_", 1:4), paste0("8cell_", 1:8))
mat <- list()
for (i in 1:length(stage)) {
  poster <- read.csv(paste0(stage[i], "_enhancer_tss2kb_posterior.txt"), sep = "\t")
  mat[[i]] <- apply(poster,1,mean)
}
all_mat <- do.call(rbind,mat)
write.table(mat, "all_enhancer_tss2kb_posterior.txt", sep = "\t", quote = F)

############################ preprocessing #########################
library(matrixStats)
mat <- read.csv("all_enhancer_tss2kb_posterior.txt", sep = "\t" )
mat <- mat[,-1]
mat <- mat[which(rowSums(mat) > 0),]
mat_rou <- round(mat, 2)
rownames(mat_rou) <- rownames(mat)
mat_rou <- mat_rou[which(rowSums(mat_rou) > 0),]
mat_rou <- as.matrix(mat_rou)
mat_rou <- mat_rou[which(rowSds(mat_rou) >0),]
mat_rou <- as.data.frame(mat_rou)

# combine rows that are exactly the same for all columns
is_equal_to_previous <- c(FALSE, apply(mat_rou[-1, ] == mat_rou[-nrow(mat_rou), ], 1, all))
deleted_rows <- numeric(length = nrow(mat_rou))
start_range <- 1
for (i in 2:nrow(mat_rou)) {
  if (!is_equal_to_previous[i]) {
    deleted_rows[start_range] <- sum(is_equal_to_previous[start_range:(i - 1)])
    start_range <- i
  }
}
filtered_mat_rou <- mat_rou[!is_equal_to_previous, ]
filtered_mat_rou$deleted_rows <- deleted_rows[!is_equal_to_previous]
rownames(filtered_mat_rou) <- rownames(mat_rou)[!is_equal_to_previous]

filtered_mat_rou <- filtered_mat_rou[filtered_mat_rou$deleted_rows != 0,]
filtered_mat_rou <- filtered_mat_rou[filtered_mat_rou$deleted_rows != 1,]

write.table(filtered_mat_rou, "all_enhancer_tss2kb_posterior_selected.txt", sep = "\t", quote = F)





