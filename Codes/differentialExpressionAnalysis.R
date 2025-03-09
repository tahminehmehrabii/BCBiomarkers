library(limma)
library(readr)
library(data.table)

# Set the working directory
setwd("C:/Users/MHR/Desktop/PDAC")

# Log2 transform
log2trans <- function(expr) {
  quan <- quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE)
  if ((quan[5] > 100) || (quan[6] - quan[1] > 50 && quan[2] > 0)) {
    expr[expr <= 0] <- NaN
    expr <- log2(expr)
  }
  expr
}

# Differential expression analysis
degs <- function(expr, group) {
  design.mat <- model.matrix(~ group + 0)
  colnames(design.mat) <- levels(group)
  fit <- lmFit(expr, design.mat)
  cont.mat <- makeContrasts(Tumor-Normal, levels = design.mat)
  fit <- eBayes(contrasts.fit(fit, cont.mat), 0.01)
  topTable(fit, number = Inf, adjust.method = "fdr", sort.by = "B")
}

# Find up & down regulated genes
updown <- function(tT) {
  up.reg <- subset(tT, adj.P.Val < 0.05 & logFC > 1)
  down.reg <- subset(tT, adj.P.Val < 0.05 & logFC < -1)
  write_tsv(up.reg, "topTable.up.tsv")
  write_tsv(down.reg, "topTable.down.tsv")
  write.table(unique(up.reg$Gene.symbol), "up.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(unique(down.reg$Gene.symbol), "down.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(union(up.reg$Gene.symbol, down.reg$Gene.symbol), "updown.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Load dataset and process
trainingSet <- fread("corrected_training_set.csv")
group <- factor(gsub(c("1" = "N", "2" = "C"), as.character(trainingSet$group), levels = c("Tumor", "Normal")))
rownames(trainingSet) <- trainingSet$V1
trainingSet <- as.data.frame(t(trainingSet[, -c(1, 2)]))

# Log2 transform and perform DE analysis
trainingSet <- log2trans(trainingSet)
tT <- degs(trainingSet, group)

# Format and save results
tT$Gene.symbol <- rownames(tT)
tT <- tT[, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
write_tsv(tT, "topTable.tsv")
write_csv(tT, "topTable.csv")

# Save up/down regulated genes
updown(tT)
