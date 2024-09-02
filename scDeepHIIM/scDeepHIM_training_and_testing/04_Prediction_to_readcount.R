Sample1_prediction<-read.table("")  #Prediction file (files needs to be transposed)
Sample1_true<-read.table("")   #Ground truth file (files needs to be transposed)
Sample1<-cbind(Sample1_prediction,Sample1_true)
colnames(Sample1)<-c("Predicted_H3K4me3","Predicted_H3K4me1","Predicted_H3K36me3","True_H3K4me3","True_H3K4me1","True_H3K36me3")

library(ggcorrplot)
# plot corr heatmap
pdf("Correlation_Heatmap.pdf")
correlation<-cor(Sample1,method = "pearson")
ggcorrplot::ggcorrplot(correlation,hc.order = T)
dev.off()

Matrix_10kb<-read.table("mm10.10K.windows.txt")
Prediction_readCount_matrix_H3K36me3<-cbind(Matrix_10kb[,c(1:3)],Sample1_prediction[,1])
Prediction_readCount_matrix_H3K4me1<-cbind(Matrix_10kb[,c(1:3)],Sample2_prediction[,2])
Prediction_readCount_matrix_H3K4me3<-cbind(Matrix_10kb[,c(1:3)],Sample3_prediction[,3])
Prediction_readCount_matrix_H3K36me3[Prediction_readCount_matrix_H3K36me3<0] <- 0
Prediction_readCount_matrix_H3K4me1[Prediction_readCount_matrix_H3K4me1<0] <- 0
Prediction_readCount_matrix_H3K4me3[Prediction_readCount_matrix_H3K4me3<0] <- 0

index1<-10000000/sum(Sample1_prediction[,1])
index2<-10000000/sum(Sample1_prediction[,2])
index3<-10000000/sum(Sample1_prediction[,3])
index4<-10000000/sum(Sample1_true[,1])
index5<-10000000/sum(Sample1_true[,2])
index6<-10000000/sum(Sample1_true[,3])

True_readCount_matrix_H3K36me3<-cbind(Matrix_10kb[,c(1:3)],Sample1_true[,1])
True_readCount_matrix_H3K4me1<-cbind(Matrix_10kb[,c(1:3)],Sample1_true[,2])
True_readCount_matrix_H3K4me3<-cbind(Matrix_10kb[,c(1:3)],Sample1_true[,3])

Prediction_readCount_matrix_H3K36me3[,4]<-Prediction_readCount_matrix_H3K36me3[,4]*index1
Prediction_readCount_matrix_H3K4me1[,4]<-Prediction_readCount_matrix_H3K4me1[,4]*index2
Prediction_readCount_matrix_H3K4me3[,4]<-Prediction_readCount_matrix_H3K4me3[,4]*index3
True_readCount_matrix_H3K36me3[,4]<-True_readCount_matrix_H3K36me3[,4]*index4
True_readCount_matrix_H3K4me1[,4]<-True_readCount_matrix_H3K4me1[,4]*index5
True_readCount_matrix_H3K4me3[,4]<-True_readCount_matrix_H3K4me3[,4]*index6

Prediction_readCount_matrix_H3K36me3[,4]<-round(Prediction_readCount_matrix_H3K36me3[,4])
Prediction_readCount_matrix_H3K4me1[,4]<-round(Prediction_readCount_matrix_H3K4me1[,4])
Prediction_readCount_matrix_H3K4me3[,4]<-round(Prediction_readCount_matrix_H3K4me3[,4])
True_readCount_matrix_H3K36me3[,4]<-round(True_readCount_matrix_H3K36me3[,4])
True_readCount_matrix_H3K4me1[,4]<-round(True_readCount_matrix_H3K4me1[,4])
True_readCount_matrix_H3K4me3[,4]<-round(True_readCount_matrix_H3K4me3[,4])

write.table(Prediction_readCount_matrix_H3K36me3, file = "Prediction_readCount_matrix_H3K36me3.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(Prediction_readCount_matrix_H3K4me1, file = "Prediction_readCount_matrix_H3K4me1.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(Prediction_readCount_matrix_H3K4me3, file = "Prediction_readCount_matrix_H3K4me3.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(True_readCount_matrix_H3K36me3, file = "True_readCount_matrix_H3K36me3.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(True_readCount_matrix_H3K4me1, file = "True_readCount_matrix_H3K4me1.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(True_readCount_matrix_H3K4me3, file = "True_readCount_matrix_H3K4me3.txt", sep = "\t", quote = FALSE, row.names = FALSE,col.names = FALSE)
