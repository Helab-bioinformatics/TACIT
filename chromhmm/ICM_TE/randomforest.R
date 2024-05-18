
##########################################################################################
library(tidyr)
library(matrixStats)

all_mat <- read.csv("./all_cells_states.txt", sep = " ")
all_mat <- all_mat[,grep("ICM|TE", colnames(all_mat))]
all_mat[all_mat=="pr"] <- 1
all_mat[all_mat=="en"] <- 2
all_mat[all_mat=="ge"] <- 3
all_mat[all_mat=="he"] <- 4
all_mat[all_mat=="qu"] <- 5

all_mat<- apply(all_mat, 2, as.numeric)
rownames(all_mat) <- rownames(all_dat_2)
all_mat <- as.matrix(all_mat)
all_mat <- all_mat[which(rowSds(all_mat) >0),]
all_mat <- all_mat[which(rowSums(all_mat == 0) <11),]

################################# randomForest ###############################
library(randomForest)
library(pROC)
library(ROCR)
library(matrixStats)
library(tidyr)
set.seed(12)

all_mat <- as.data.frame(t(all_mat))
all_mat$identity <- as.factor(c(rep("ICM", 15), rep("TE", 24)))
index <- sample(1:39, 3)
train <- all_mat[-index,]
test <- all_mat[index,]
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

pdf("./train.pdf")
plot(mtry_fit)#best ntree=10000
plot(mtry_fit$err.rate)
title(main = "use all diff bins", sub = "err=11.11%")


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

#filter importance_otu
importance_otu <- importance(mtry_fit)
varImpPlot(mtry_fit)
importance_otu <- as.data.frame(importance_otu)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(importance_otu)
#write.table(importance_otu, './randomForest/importance_otu_top1500.txt', sep = '\t', col.names = NA, quote = FALSE)
#write.table(importance_otu[1:30, ], 'importance_otu_top30.txt', sep = '\t', col.names = NA, quote = FALSE)

set.seed(123)
otu_train.cv <- replicate(5, rfcv(train[-ncol(train)], train$identity, cv.fold = 10,step = 1.5), simplify = FALSE)
otu_train.cv

otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))

library(ggplot2)
library(splines) 

p <- ggplot(otu_train.cv, aes(otus, value)) +
  geom_smooth(se = FALSE,  method = 'glm', formula = y~ns(x, 6)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  labs(title = '',x = 'Number of OTUs', y = 'Cross-validation error')
p
dev.off()

