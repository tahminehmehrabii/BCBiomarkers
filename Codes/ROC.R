library(pROC)
library(data.table)

set.seed(123)

# Set the current working directory to the project path
setwd("C:/Users/MHR/Desktop/Breast")

#### ROC curve function ####
rocCurve <- function(data, genes, path, mfrow, auc.x, auc.y, w, h) {
  png(filename = path, width = w, height = h, res = 300)
  par(mfrow = mfrow)
  for(g in genes) {
    roc.train <- roc(as.formula(paste("group ~", g)), data = data[[1]], auc = TRUE, ci = TRUE)
    roc.valid <- roc(as.formula(paste("group ~", g)), data = data[[2]], auc = TRUE, ci = TRUE)
    
    #ROC
    plot.roc(roc.train, 
             col = "#6A0DAD",  
             legacy.axes = FALSE,
             main = g,
             cex.main = 1.3,
             lwd = 2.3,
             cex.lab = 1.4,
             type = "s"   
    )
    
    lines(roc.valid, 
          col = "#FF1493",   
          lwd = 2.3,
          type = "s"         
    )
    
    legend("bottomright",
           legend = c(paste0("train AUC = ", signif(roc.train$auc * 100, digits = 4), "%"),
                      paste0("valid AUC = ", signif(roc.valid$auc * 100, digits = 4), "%")),
           col = c("#6A0DAD", "#FF1493"),
           lty = 1,
           lwd = 2.2,
           cex = 1.5
    )
  }
  dev.off()
}


# Load selected biomarkers
biomarkers <- fread("lassoGenes.txt", header = FALSE)$V1
biomarkers <- sort(biomarkers)

#### Load training data ####
trainingSet <- as.data.frame(fread("corrected_training_set.csv"))
rownames(trainingSet) <- trainingSet$sample
group <- as.factor(trainingSet$group)
levels(group) <- c("Normal", "Tumor")
trainingSet <- trainingSet[ , -c(1:4)]
trainingSet <- trainingSet[ , biomarkers, drop = FALSE]
train <- cbind(group = group, trainingSet)

#### Load validation data ####
validSet <- as.data.frame(fread("corrected_validation_set.csv"))
rownames(validSet) <- validSet$sample
group <- as.factor(validSet$group)
levels(group) <- c("Normal", "Tumor")
validSet <- validSet[ , -c(1:4)]
validSet <- validSet[ , biomarkers, drop = FALSE]
valid <- cbind(group = group, validSet)

# Combine data
data <- list(train, valid)

#### Filter genes with AUC > 0.92 in train and > 0.92 in valid sets ####
selected_genes <- biomarkers[sapply(biomarkers, function(g) {
  roc.train <- roc(as.formula(paste("group ~", g)), data = data[[1]], auc = TRUE, ci = TRUE)
  roc.valid <- roc(as.formula(paste("group ~", g)), data = data[[2]], auc = TRUE, ci = TRUE)
  (roc.train$auc > 0.92) & (roc.valid$auc > 0.92)
})]


# Auto-calculate mfrow dimensions based on number of genes
rows <- ceiling(sqrt(length(selected_genes)))
mfrow_dims <- c(rows, rows)

# Run ROC curve function for selected genes
rocCurve(data, selected_genes, "Filtered_ROC.png", mfrow_dims, 0.78, 0.50, 3400, 3400)

# Save processed datasets to CSV
fwrite(train, "training_data.csv", row.names = FALSE)
fwrite(valid, "validation_data.csv", row.names = FALSE)
