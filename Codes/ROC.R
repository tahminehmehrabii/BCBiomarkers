library(pROC)
library(data.table)

set.seed(123)

# Set the current working directory to the project path
setwd("C:/Users/MHR/Desktop/Breast")

#### ROC curve function ####
rocCurve <- function(data, genes, path, mfrow, auc.x, auc.y, w, h) {
  png(filename = path, width = w, height = h, res = 300)
  par(mfrow=mfrow)
  for(g in genes) {
    roc.train <- roc(as.formula(paste("group ~", g)), data = data[[1]], auc = TRUE, ci = TRUE)
    roc.valid <- roc(as.formula(paste("group ~", g)), data = data[[2]], auc = TRUE, ci = TRUE)
    plot.roc(roc.train, 
             col = "#1F77B4",
             legacy.axes = FALSE,
             main = g,
             cex.main = 1.3,
             lwd = 2.3,
             cex.lab = 1.4
    )
    lines(roc.valid, col="salmon")
    legend("bottomright",
           legend=c(paste0("train AUC=", signif(roc.train$auc*100, digits = 4), "%"),
                    paste0("valid AUC=", signif(roc.valid$auc*100, digits = 4), "%")),
           col=c("#1F77B4", "salmon"),
           lty = 1,
           lwd = 2.2,
           cex = 1.5
    )
  }
  dev.off()
}

biomarkers <- fread("lassoGenes.txt", header = FALSE)$V1
biomarkers <- sort(biomarkers)

#### ROC curve for train data ####
trainingSet <- as.data.frame(fread("corrected_training_set.csv"))
rownames(trainingSet) <- trainingSet$sample
group <- trainingSet$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
trainingSet <- trainingSet[,-c(1:4)]
trainingSet <- trainingSet[, biomarkers]
train <- cbind(group=group, trainingSet)

#### ROC curve for testing data ####
validSet <- as.data.frame(fread("corrected_validation_set.csv"))
rownames(validSet) <- validSet$sample
group <- validSet$group
group <- as.factor(group)
levels(group) <- c("Normal", "Tumor")
validSet <- validSet[,-c(1:4)]
validSet <- validSet[, biomarkers]
valid <- cbind(group=group, validSet)

data <- list(train, valid)

# Filter genes with AUC > 0.9 in both sets
selected_genes <- biomarkers[sapply(biomarkers, function(g) {
  roc.train <- roc(as.formula(paste("group ~", g)), data = data[[1]], auc = TRUE, ci = TRUE)
  roc.valid <- roc(as.formula(paste("group ~", g)), data = data[[2]], auc = TRUE, ci = TRUE)
  (roc.train$auc > 0.9) & (roc.valid$auc > 0.9)
})]

# Run ROC curve function for selected genes
rocCurve(data, selected_genes, "Filtered_ROC.png", c(3, 3), 0.78, 0.50, 3400, 3400)



# Save processed datasets to CSV
fwrite(train, "training_data.csv", row.names = FALSE)
fwrite(valid, "validation_data.csv", row.names = FALSE)
