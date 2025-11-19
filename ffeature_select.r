library(caret)
library(glmnet)
library(doParallel)
library(foreach)
library(dplyr)
library(pROC)
#setwd("D:/02.single_cell_screen/06machine_learning")
df <- read.table("log2fc_of_treat_vs_CK_pairs_combined_135.20250102.xls",sep = "\t",stringsAsFactors = F,header = T,row.names = 136)
df <- as.data.frame(t(df))
winfo <- read.table("info_of_treat_vs_CK_pairs_combined_135.20250102.xls",sep = "\t",stringsAsFactors = F,header = T)
winfo$X3 <- tolower(winfo$X3)
rownames(winfo) <- winfo$name
df$label <- winfo[rownames(df),"X3"]
df$label <- factor(df$label)
registerDoParallel(cores = 30)
results <- foreach(i=c(695),.packages = c("glmnet","caret","dplyr","pROC")) %dopar% {
  set.seed(i)
  trainIndex <- createDataPartition(df$label, p = 0.8, list = FALSE)
  trainData0 <- df[trainIndex,-ncol(df)]
  trainData0_var <- apply(trainData0,2,var)
  trainData0 <- trainData0[,trainData0_var>0.01]
  trainData0_cor <- cor(trainData0)
  highly_correlated <- findCorrelation(trainData0_cor, cutoff = 0.75)
  trainData0 <- trainData0[, -highly_correlated]
  testData0 <- df[-trainIndex,colnames(trainData0)]
  preProc <- preProcess(trainData0, method = c("center", "scale"))
  trainData <- predict(preProc, trainData0)
  testData <- predict(preProc, testData0)
  ###lasso 特征选择
  trainData_lasso <- trainData
  trainData_lasso$label <- winfo[rownames(trainData_lasso),"X3"]
  trainData_lasso$label <- ifelse(trainData_lasso$label=="tolerant",1,0)
  x <- model.matrix(label ~ ., trainData_lasso)[,-1]
  y <- trainData_lasso$label
  set.seed(123)
  cv_lasso <- cv.glmnet(x, y, alpha = 1,nfolds = 5)
  best_lambda <- cv_lasso$lambda.min
  lasso_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
  coef <- coef(lasso_model)
  selected_features <- rownames(coef)[coef[, 1] != 0][-1]
  if(length(selected_features)>=2){
    trainData2 <- trainData0[,selected_features]
    trainData2$label <- factor(winfo[rownames(trainData2),"X3"])
    testData2 <- testData0[,selected_features]
    testData2$label <- factor(winfo[rownames(testData2),"X3"])
    
    preProc <- preProcess(trainData2, method = c("center", "scale"))
    trainData <- predict(preProc, trainData2)
    testData <- predict(preProc, testData2)
    save(preProc,file = paste("seed",i,"_scale.RData",sep = ""))
    
    train.control <- trainControl(method = "repeatedcv",
                                  number = 5,
                                  repeats = 10,
                                  classProbs = TRUE,
                                  savePredictions = "final",
                                  summaryFunction=twoClassSummary)
    glmnetfit <- train(x = trainData[,-ncol(trainData)],
                       y= trainData$label,
                       method = "ranger",
                       trControl = train.control,
                       tuneLength = 16,
                       metric = "ROC")
    fit <- glmnetfit
    test_pred <- select(testData,label) %>%
      bind_cols(predict(fit,newdata=testData,type="prob")) %>%
      bind_cols(pred=predict(fit,newdata=testData))
    roc <- roc(response = test_pred$label, predictor = test_pred$tolerant, levels = rev(levels(test_pred$label)))
    cfm <- confusionMatrix(data=test_pred$label,test_pred$pred,positive = "tolerant")
    auc <- as.numeric(roc$auc)
    accuracy <- as.numeric(cfm$overall["Accuracy"])
  }else{
    auc <- NA
    accuracy <- NA
  }
  r <- c(i,length(selected_features),min(cv_lasso$cvm),auc,accuracy,paste(selected_features,collapse = ","))
  r <- paste(r,collapse = ";")
  write.table(r,paste("seed",i,"_feature_MSE.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  coef_df <- data.frame(id=coef@Dimnames[[1]],coef=coef[,1])
  coef_df <- coef_df[-1,]
  coef_df <- coef_df[coef_df$coef!=0,]
  write.table(coef_df,paste("seed",i,"_features.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
  save(fit,file=paste("seed",i,"_ranger_model.RData",sep = ""))
}

stopImplicitCluster()
