library(caret)
library(glmnet)
library(doParallel)
library(foreach)
library(dplyr)
library(pROC)
#setwd("D:/02.single_cell_screen/06machine_learning/01rice_wheat_lg2fc_135_samples/")
df <- read.table("log2fc_of_treat_vs_CK_pairs_combined_135.20250102.xls",sep = "\t",stringsAsFactors = F,header = T,row.names = 136)
df <- as.data.frame(t(df))
winfo <- read.table("info_of_treat_vs_CK_pairs_combined_135.20250102.xls",sep = "\t",stringsAsFactors = F,header = T)
winfo$X3 <- tolower(winfo$X3)
rownames(winfo) <- winfo$name
df$label <- winfo[rownames(df),"X3"]
df$label <- factor(df$label)
registerDoParallel(cores = 15)
results <- foreach(j=c(731),.packages = c("glmnet","caret","dplyr","pROC")) %dopar% {
  i <- 411
  set.seed(i)
  trainIndex <- createDataPartition(df$label, p = 0.8, list = FALSE)
  select_features <- read.table(paste("glmnet/seed",i,"_features.txt",sep = ""),stringsAsFactors = F,sep = "\t",header = T)
  select_features <- select_features$id
  trainData0 <- df[trainIndex,select_features]
  testData0 <- df[-trainIndex,select_features]
  #load(paste("glmnet/seed",i,"_scale.RData",sep = ""))
  #load(paste("seed",i,"_ranger_model.RData",sep = ""))
  #preProc <- preProcess(trainData0, method = c("center", "scale"))
  #trainData <- predict(preProc, trainData0)
  #testData <- predict(preProc, testData0)
  trainData <- trainData0
  testData <- testData0
  trainData$label <- factor(winfo[rownames(trainData),"X3"],levels=c("tolerant","sensitive"))
  testData$label <- factor(winfo[rownames(testData),"X3"],levels=c("tolerant","sensitive"))
  train.control <- trainControl(method = "repeatedcv",
                                number = 5,
                                repeats = 10,
                                classProbs = TRUE,
                                savePredictions = "final",
                                summaryFunction=twoClassSummary)
  set.seed(j)
  fit <- train(x = trainData[,-ncol(trainData)],
               y= trainData$label,
               method = "ranger",
               trControl = train.control,
               tuneLength = 16,
	       importance = "impurity",
               metric = "ROC")
  save(fit,file=paste("seed411_seed",j,"_ranger_modelv3.RData",sep = ""))
  test_pred <- select(testData,label) %>%
    bind_cols(predict(fit,newdata=testData[,-ncol(testData)],type="prob")) %>%
    bind_cols(pred=predict(fit,newdata=testData[,-ncol(testData)]))
  roc <- roc(response = test_pred$label, predictor = test_pred$tolerant, levels = rev(levels(test_pred$label)))
  cfm <- confusionMatrix(data=test_pred$label,test_pred$pred,positive = "tolerant")
  auc <- as.numeric(roc$auc)
  accuracy <- as.numeric(cfm$overall["Accuracy"])
  preci <- as.numeric(cfm$byClass["Precision"])
  F1 <- as.numeric(cfm$byClass["F1"])
  recal <- as.numeric(cfm$byClass["Recall"])
  
  r <- c(411,j,auc,accuracy,preci,F1,recal)
  r <- paste(r,collapse = ";")
  write.table(r,paste("seed411_seed",j,"_model_preformancev3.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
}

