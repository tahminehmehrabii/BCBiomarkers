# GSE29431*********************************************
#########################################################
library(data.table)
library(ggplot2)
library(caret)
library(pROC)
library(MLmetrics)
library(randomForest)
library(e1071)
library(neuralnet)
library(nnet)
library(dplyr)
library(NeuralNetTools)
library(h2o)
library(GEOquery)

setwd("C:/Users/MHR/Desktop/Breast")

set.seed(123)

get_exprs_matrix <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE)[[1]]
  return(list(expr = exprs(gse), gse = gse))
}

data29431 <- get_exprs_matrix("GSE29431")
expr29431 <- data29431$expr
gse29431 <- data29431$gse


write.table(expr29431, file = "GSE29431.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

feat29431 <- fData(gse29431)
write.table(feat29431, file = "feat29431.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

phen29431 <- pData(gse29431)
write.table(phen29431, file = "phen29431.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

expr29431 <- cbind(ID=rownames(expr29431), expr29431)
write.table(expr29431, file = "expr29431.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

feat29431 <- cbind(ID = feat29431[, "ID"], symbol = feat29431[, "Gene Symbol"])

expr29431<- merge(feat29431, expr29431, by="ID")

rownames(expr29431) <- expr29431$ID

expr29431 <- expr29431[,-1]

lassoGenes <- fread("lassoGenes.txt", header = FALSE)$V1
lassoGenes <- sort(lassoGenes)

biomarkers <- lassoGenes

expr29431 <- expr29431[which(expr29431$symbol %in% biomarkers),]

rownames(expr29431) <- NULL

expr29431 <- expr29431[!duplicated(expr29431$symbol),]


rownames(expr29431) <- expr29431$symbol

expr29431 <- expr29431[,-1]

expr29431 <- as.data.frame(t(expr29431))

expr29431 <- expr29431[,biomarkers]


for (j in 1:ncol(expr29431)) expr29431[,j] <- as.numeric(expr29431[,j])

expr29431 <- na.omit(expr29431)

expr29431 <- cbind(sample=rownames(expr29431), expr29431)

phen29431 <- cbind(sample=rownames(phen29431), phen29431)


phen29431  <- data.frame(
  sample = phen29431 [, "geo_accession"], 
  group = phen29431 [, "characteristics_ch1"]  
)
phen29431 $group <- ifelse(phen29431 $group == "disease state: none", "Normal", "Tumor")




expr29431<- merge(phen29431, expr29431, by="sample")

rownames(expr29431) <- expr29431$sample

expr29431 <- expr29431[,-1]

write.csv(expr29431, "trainEvaluateML_GSE29431.csv")
#############################################################################################
##############################################
trainSet <- as.data.frame(fread("corrected_training_set.csv"))
validSet <- as.data.frame(fread("corrected_validation_set.csv"))
gse29431 <- as.data.frame(fread("trainEvaluateML_GSE29431.csv"))

trainGroup <- trainSet$group
validGroup <- validSet$group
gse29431Group <- gse29431$group

rownames(trainSet) <- trainSet$sample
rownames(validSet) <- validSet$sample
rownames(gse29431) <- gse29431$V1

trainData <- trainSet[,-c(1:4)]
validData <- validSet[,-c(1:4)]
gse29431Data <- gse29431[,-c(1,2)]

biomarkers <- fread("lassoGenes.txt", header = FALSE)$V1
biomarkers <- sort(biomarkers)

trainData <- trainData[,which(colnames(trainData) %in% biomarkers)]
validData <- validData[,which(colnames(validData) %in% biomarkers)]
gse29431Data <- gse29431Data[,which(colnames(gse29431Data) %in% biomarkers)]

trainData <- trainData[,biomarkers]
validData <- validData[,biomarkers]
gse29431Data <- gse29431Data[,biomarkers]

trainGroup <- as.factor(trainGroup)
validGroup <- as.factor(validGroup)
gse29431Group <- as.factor(gse29431Group)

levels(trainGroup) <- c("Normal", "Tumor")
levels(validGroup) <- c("Normal", "Tumor")
levels(gse29431Group) <- c("Normal", "Tumor")

train <- cbind(group=trainGroup, trainData)
valid <- cbind(group=validGroup, validData)
gse29431 <- cbind(group=gse29431Group, gse29431Data)

test <- rbind(gse29431)
accession <- "GSE29431"

write.csv(train, "train_data_GSE29431.csv")
write.csv(valid, "valid_data_GSE29431.csv")
write.csv(test, "test_data_GSE29431.csv")




###################################

##############################################

library(randomForest)
library(pROC)
library(caret)
library(MLmetrics)

#### Random forest ####
#######################

ntree = 500

tuned_rf <- tuneRF(train[,-1], 
                   train$group,
                   ntreeTry = ntree,
                   stepFactor = 1.5,
                   improve = 0.05,
                   trace = TRUE,
                   plot = FALSE
)

best.m <- tuned_rf[tuned_rf[,2] == min(tuned_rf[,2]),1][1]

best.rf <- randomForest(formula = group ~ ., 
                        data = train, 
                        mtry=best.m,
                        importance=TRUE,
                        ntree=ntree
)

rf.train.group.true <- train$group
rf.train.group.pred <- as.factor(predict(best.rf, train[,-1]))
rf.train.group.prob <- as.data.frame(predict(best.rf, train[,-1], type = "prob"))

rf.train.roc <- roc(rf.train.group.true ~ as.numeric(rf.train.group.prob$Tumor), ci = TRUE)

rf.train.cm <- confusionMatrix(rf.train.group.true, rf.train.group.pred, positive = "Tumor")
write.csv(rf.train.cm$table, "rf.train.confusion.matrix_GSE29431.csv")

rf.train.group.true <- as.numeric(rf.train.group.true)
rf.train.group.pred <- as.numeric(rf.train.group.pred)

rf.train.matrics <- data.frame("AUC"=rf.train.roc$auc[1],
                               "CI"=paste0(round(rf.train.roc$ci[1], 4), "-", round(rf.train.roc$ci[3], 4)),
                               "accuracy"=Accuracy(rf.train.group.pred, rf.train.group.true),
                               "sensitivity"=Sensitivity(rf.train.group.true, rf.train.group.pred, positive = 2),
                               "specificity"=Specificity(rf.train.group.true, rf.train.group.pred, positive = 2),
                               "precision"=Precision(rf.train.group.true, rf.train.group.pred, positive = 2),
                               "recall"=Recall(rf.train.group.true, rf.train.group.pred, positive = 2),
                               "f1score"=F1_Score(rf.train.group.true, rf.train.group.pred, positive = 2)
)

write.csv(rf.train.matrics, "rf.train.matrics_GSE29431.csv")

rf.valid.group.true <- valid$group
rf.valid.group.pred <- as.factor(predict(best.rf, valid[,-1]))
rf.valid.group.prob <- as.data.frame(predict(best.rf, valid[,-1], type = "prob"))

rf.valid.roc <- roc(rf.valid.group.true ~ as.numeric(rf.valid.group.prob$Tumor), ci = TRUE)

rf.valid.cm <- confusionMatrix(rf.valid.group.true, rf.valid.group.pred, positive = "Tumor")
write.csv(rf.valid.cm$table, "rf.valid.confusion.matrix_GSE29431.csv")

rf.valid.group.true <- as.numeric(rf.valid.group.true)
rf.valid.group.pred <- as.numeric(rf.valid.group.pred)

rf.valid.matrics <- data.frame("AUC"=rf.valid.roc$auc[1],
                               "CI"=paste0(round(rf.valid.roc$ci[1], 4), "-", round(rf.valid.roc$ci[3], 4)),
                               "accuracy"=Accuracy(rf.valid.group.pred, rf.valid.group.true),
                               "sensitivity"=Sensitivity(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                               "specificity"=Specificity(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                               "precision"=Precision(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                               "recall"=Recall(rf.valid.group.true, rf.valid.group.pred, positive = 2),
                               "f1score"=F1_Score(rf.valid.group.true, rf.valid.group.pred, positive = 2)
)

write.csv(rf.valid.matrics, "rf.valid.matrics_GSE29431.csv")

rf.test.group.true <- test$group
rf.test.group.pred <- as.factor(predict(best.rf, test[,-1]))
rf.test.group.prob <- as.data.frame(predict(best.rf, test[,-1], type = "prob"))

rf.test.roc <- roc(rf.test.group.true ~ as.numeric(rf.test.group.prob$Tumor), ci = TRUE)

rf.test.cm <- confusionMatrix(rf.test.group.true, rf.test.group.pred, positive = "Tumor")
write.csv(rf.test.cm$table, paste0("rf.test.confusion.matrix_GSE29431.csv"))

rf.test.group.true <- as.numeric(rf.test.group.true)
rf.test.group.pred <- as.numeric(rf.test.group.pred)

rf.test.matrics <- data.frame("AUC"=rf.test.roc$auc[1],
                              "CI"=paste0(round(rf.test.roc$ci[1], 4), "-", round(rf.test.roc$ci[3], 4)),
                              "accuracy"=Accuracy(rf.test.group.pred, rf.test.group.true),
                              "sensitivity"=Sensitivity(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "specificity"=Specificity(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "precision"=Precision(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "recall"=Recall(rf.test.group.true, rf.test.group.pred, positive = 2),
                              "f1score"=F1_Score(rf.test.group.true, rf.test.group.pred, positive = 2)
)

write.csv(rf.test.matrics, paste0("rf.test.matrics_GSE29431.csv"))
#######################
library(pROC)
library(caret)
library(MLmetrics)
library(e1071)

#### Support Vector Machine ####
################################

tune.svm <- tune(svm, group ~ ., 
                 data = train, 
                 kernel = "linear",
                 ranges = list(cost = c(0.1, 1, 10),
                               gamma = c(0.1, 0.01, 0.001))
)

svm.fit <- svm(group ~ .,
               data = train,
               kernel = "linear",
               cost = tune.svm$best.parameters$cost,
               gamma = tune.svm$best.parameters$gamma,
               probability = TRUE
)

# Prediction on train data
svm.train.group.true <- train$group
svm.train.group.pred <- predict(svm.fit, newdata = train[,-1])
svm.train.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = train[,-1], probability = TRUE), "probabilities"))
svm.train.roc <- roc(svm.train.group.true, svm.train.group.prob$Tumor, ci = TRUE)
svm.train.cm <- confusionMatrix(svm.train.group.true, svm.train.group.pred, positive = "Tumor")
write.csv(svm.train.cm$table, "svm.train.confusion.matrix_GSE29431.csv")

svm.train.matrics <- data.frame("AUC"=svm.train.roc$auc[1],
                                "CI"=paste0(round(svm.train.roc$ci[1], 4), "-", round(svm.train.roc$ci[3], 4)),
                                "accuracy"=Accuracy(svm.train.group.pred, svm.train.group.true),
                                "sensitivity"=Sensitivity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                                "specificity"=Specificity(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                                "precision"=Precision(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                                "recall"=Recall(svm.train.group.true, svm.train.group.pred, positive = "Tumor"),
                                "f1score"=F1_Score(svm.train.group.true, svm.train.group.pred, positive = "Tumor")
)

write.csv(svm.train.matrics, "svm.train.matrics_GSE29431.csv")

# Prediction on validation data
svm.valid.group.true <- valid$group
svm.valid.group.pred <- predict(svm.fit, newdata = valid[,-1])
svm.valid.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = valid[,-1], probability = TRUE), "probabilities"))
svm.valid.roc <- roc(svm.valid.group.true, svm.valid.group.prob$Tumor, ci = TRUE)
svm.valid.cm <- confusionMatrix(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor")
write.csv(svm.valid.cm$table, "svm.valid.confusion.matrix_GSE29431.csv")

svm.valid.matrics <- data.frame("AUC"=svm.valid.roc$auc[1],
                                "CI"=paste0(round(svm.valid.roc$ci[1], 4), "-", round(svm.valid.roc$ci[3], 4)),
                                "accuracy"=Accuracy(svm.valid.group.pred, svm.valid.group.true),
                                "sensitivity"=Sensitivity(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                                "specificity"=Specificity(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                                "precision"=Precision(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                                "recall"=Recall(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor"),
                                "f1score"=F1_Score(svm.valid.group.true, svm.valid.group.pred, positive = "Tumor")
)

write.csv(svm.valid.matrics, "svm.valid.matrics_GSE29431.csv")

# Prediction on test data
svm.test.group.true <- test$group
svm.test.group.pred <- predict(svm.fit, newdata = test[,-1])
svm.test.group.prob <- as.data.frame(attr(predict(svm.fit, newdata = test[,-1], probability = TRUE), "probabilities"))
svm.test.roc <- roc(svm.test.group.true, svm.test.group.prob$Tumor, ci = TRUE)
svm.test.cm <- confusionMatrix(svm.test.group.true, svm.test.group.pred, positive = "Tumor")
write.csv(svm.test.cm$table, "svm.test.confusion.matrix_GSE29431.csv")

svm.test.matrics <- data.frame("AUC"=svm.test.roc$auc[1],
                               "CI"=paste0(round(svm.test.roc$ci[1], 4), "-", round(svm.test.roc$ci[3], 4)),
                               "accuracy"=Accuracy(svm.test.group.pred, svm.test.group.true),
                               "sensitivity"=Sensitivity(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "precision"=Precision(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "recall"=Recall(svm.test.group.true, svm.test.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(svm.test.group.true, svm.test.group.pred, positive = "Tumor")
)

write.csv(svm.test.matrics, "svm.test.matrics_GSE29431.csv")


##############################################################################################
#### Artificial Neural Network ####
###################################

# Encode as a one hot vector multilabel data
train.onehot <- cbind(class.ind(as.factor(train$group)), train[,-1])
valid.onehot <- cbind(class.ind(as.factor(valid$group)), valid[,-1])
test.onehot <- cbind(class.ind(as.factor(test$group)), test[,-1])

# Set up formula
n <- names(train.onehot)
f <- as.formula(paste("Normal + Tumor ~", paste(n[!n %in% c("Normal", "Tumor")], collapse = " + ")))

nn <- neuralnet(f,
                data = train.onehot,
                hidden = c(5),
                act.fct = "tanh",
                linear.output = FALSE)

png("annplot.png", height = 1600, width = 2200, res = 300)
plotnet(nn,
        cex_val = 0.8,
        circle_cex = 4,
        circle_col = "pink",
        bord_col = "pink",
        neg_col = "grey",
        nid = TRUE
)
dev.off()

# Macke prediction on train data
nn.train.group.prob <- as.data.frame(predict(nn, train.onehot[,3:ncol(train.onehot)]))
colnames(nn.train.group.prob) <- c("Normal", "Tumor")

nn.train.group.pred <- max.col(nn.train.group.prob) |> as.factor()
nn.train.group.true <- max.col(train.onehot[,1:2]) |> as.factor()

levels(nn.train.group.pred) <- c("Normal", "Tumor")
levels(nn.train.group.true) <- c("Normal", "Tumor")

nn.train.cm <- confusionMatrix(nn.train.group.true, nn.train.group.pred, positive = "Tumor")
write.csv(nn.train.cm$table, "nn.train.confusion.matrix_GSE29431.csv")

nn.train.roc <- roc(nn.train.group.true ~ as.numeric(nn.train.group.prob$Tumor), ci = TRUE)

nn.train.matrics <- data.frame("AUC"=nn.train.roc$auc[1],
                               "CI"=paste0(round(nn.train.roc$ci[1], 4), "-", round(nn.train.roc$ci[3], 4)),
                               "accuracy"=Accuracy(nn.train.group.pred, nn.train.group.true),
                               "sensitivity"=Sensitivity(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "precision"=Precision(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "recall"=Recall(nn.train.group.true, nn.train.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(nn.train.group.true, nn.train.group.pred, positive = "Tumor")
)

write.csv(nn.train.matrics, "nn.train.matrics_GSE29431.csv")


# Macke prediction on valid data
nn.valid.group.prob <- as.data.frame(predict(nn, valid.onehot[,3:ncol(valid.onehot)]))
colnames(nn.valid.group.prob) <- c("Normal", "Tumor")

nn.valid.group.pred <- max.col(nn.valid.group.prob) |> as.factor()
nn.valid.group.true <- max.col(valid.onehot[,1:2]) |> as.factor()

levels(nn.valid.group.pred) <- c("Normal", "Tumor")
levels(nn.valid.group.true) <- c("Normal", "Tumor")

nn.valid.cm <- confusionMatrix(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor")
write.csv(nn.valid.cm$table, "nn.valid.confusion.matrix_GSE29431.csv")

nn.valid.roc <- roc(nn.valid.group.true ~ as.numeric(nn.valid.group.prob$Tumor), ci = TRUE)

nn.valid.matrics <- data.frame("AUC"=nn.valid.roc$auc[1],
                               "CI"=paste0(round(nn.valid.roc$ci[1], 4), "-", round(nn.valid.roc$ci[3], 4)),
                               "accuracy"=Accuracy(nn.valid.group.pred, nn.valid.group.true),
                               "sensitivity"=Sensitivity(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                               "specificity"=Specificity(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                               "precision"=Precision(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                               "recall"=Recall(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor"),
                               "f1score"=F1_Score(nn.valid.group.true, nn.valid.group.pred, positive = "Tumor")
)

write.csv(nn.valid.matrics, "nn.valid.matrics_GSE29431.csv")

# Macke prediction on test data
nn.test.group.prob <- as.data.frame(predict(nn, test.onehot[,3:ncol(test.onehot)]))
colnames(nn.test.group.prob) <- c("Normal", "Tumor")

nn.test.group.pred <- max.col(nn.test.group.prob) |> as.factor()
nn.test.group.true <- max.col(test.onehot[,1:2]) |> as.factor()

levels(nn.test.group.pred) <- c("Normal", "Tumor")
levels(nn.test.group.true) <- c("Normal", "Tumor")

nn.test.cm <- confusionMatrix(nn.test.group.true, nn.test.group.pred, positive = "Tumor")
write.csv(nn.test.cm$table, paste0("nn.test.confusion.matrix_GSE29431.", accession, ".csv"))

nn.test.roc <- roc(nn.test.group.true ~ as.numeric(nn.test.group.prob$Tumor), ci = TRUE)

nn.test.matrics <- data.frame("AUC"=nn.test.roc$auc[1],
                              "CI"=paste0(round(nn.test.roc$ci[1], 4), "-", round(nn.test.roc$ci[3], 4)),
                              "accuracy"=Accuracy(nn.test.group.pred, nn.test.group.true),
                              "sensitivity"=Sensitivity(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "specificity"=Specificity(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "precision"=Precision(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "recall"=Recall(nn.test.group.true, nn.test.group.pred, positive = "Tumor"),
                              "f1score"=F1_Score(nn.test.group.true, nn.test.group.pred, positive = "Tumor")
)

write.csv(nn.test.matrics, paste0("nn.test.matrics_GSE29431.", accession, ".csv"))
######################################################################

###################################################################################
roc.curve <- function(rocs, cols, main, mods) {
  
  if (length(rocs) != length(cols) || length(rocs) != length(mods)) {
    stop("The number of models, colors, and model names must be equal!")
  }
  
  png(paste0("roc_", main ,".png"), width = 1800, height = 1800, res = 300)
  
  plot.roc(rocs[[1]],
           col = cols[1],
           legacy.axes = FALSE,
           main = main,
           cex.main = 1.5,
           lwd = 2.5,
           cex.lab = 1.5
  )
  
  for (i in 2:length(rocs)) {
    lines(rocs[[i]], col = cols[i], lwd = 2.5)
  }
  
  legend_text <- sapply(1:length(rocs), function(i) {
    paste0(mods[i], " (AUC=", signif(rocs[[i]]$auc * 100, digits = 4), "%)")
  })
  
  legend("bottomright",
         legend = legend_text,
         col = cols,
         lty = 1,
         lwd = 3,
         cex = 1.2
  )
  
  dev.off()
}

##########

###################
#### ROC plot ####
##################

cols <- c("#FF1493", "#6A0DAD", "orange")  # Changed colors

mods <- c("RF", "SVM", "ANN")

# Checking if the variables are initialized before use
if (!exists("rf.train.roc") | !exists("svm.train.roc") | !exists("nn.train.roc")) {
  stop("One of the train ROC variables is not initialized!")
}
rocs.train <- list(rf.train.roc, svm.train.roc, nn.train.roc)
roc.curve(rocs.train, cols, "Training Dataset_GSE29431", mods)

if (!exists("rf.valid.roc") | !exists("svm.valid.roc") | !exists("nn.valid.roc")) {
  stop("One of the valid ROC variables is not initialized!")
}
rocs.valid <- list(rf.valid.roc, svm.valid.roc, nn.valid.roc)
roc.curve(rocs.valid, cols, "Validation Dataset_GSE29431", mods)

# Default assignment for accession variable if not initialized
if (!exists("accession")) {
  accession <- "Unknown"
}

if (!exists("rf.test.roc") | !exists("svm.test.roc") | !exists("nn.test.roc")) {
  stop("One of the test ROC variables is not initialized!")
}
rocs.test <- list(rf.test.roc, svm.test.roc, nn.test.roc)
roc.curve(rocs.test, cols, paste0("Test Dataset (", accession, ")"), mods)

# Saving the models (checking initialization before saving)
if (exists("best.rf")) saveRDS(best.rf, "RF.rds")
if (exists("svm.fit")) saveRDS(svm.fit, "SVM.rds")
if (exists("nn")) saveRDS(nn, "ANN.rds")

