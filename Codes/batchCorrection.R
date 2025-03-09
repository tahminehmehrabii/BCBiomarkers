library(data.table)
library(RColorBrewer)
library(ggplot2)
library(sva)

# Set the current working directory to the project path
setwd("C:/Users/MHR/Desktop/PDAC")

# Load data
load_data <- function(file_path) {
  data <- as.data.frame(fread(file_path))
  rownames(data) <- data$V1
  data <- data[,-1]
  list(
    data = data[,-c(1,2)],
    batch = data$batch,
    group = data$group,
    batch_group = data[, c(1,2)]
  )
}

training <- load_data("training_set.csv")
validation <- load_data("validation_set.csv")

# Prepare variables
grps <- c("GSE71989", "GSE46234", "GSE15471")
type <- c("Normal", "Tumor")
colors <- brewer.pal(n = 5, name = "Dark2")[1:3]

# PCA function
run_pca <- function(data, batch, group, type_labels, grps_labels, colors, file_prefix) {
  pca_result <- prcomp(data, scale = FALSE, center = TRUE)
  pca_var_per <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)
  
  pc1 <- pca_result$x[,1]
  pc2 <- pca_result$x[,2]
  
  sdat <- data.frame(x = pc1, y = pc2, batch = grps_labels[batch])
  tdat <- data.frame(sdat[,1:2], grps_labels[batch], type_labels[group])
  
  colnames(tdat) <- c("PC-1", "PC-2", "Batch", "Group")
  rownames(tdat) <- rownames(data)
  write.table(tdat, paste0(file_prefix, "_pca.txt"), quote = TRUE, sep = "\t", row.names = TRUE)
  
  clrs <- c(rep(colors[1], length(which(batch == 1))),
            rep(colors[2], length(which(batch == 2))),
            rep(colors[3], length(which(batch == 3))))
  
  png(paste0(file_prefix, "_pca.png"), width = 2600, height = 2000, res = 300)
  ggplot(data = sdat, aes(x, y, colour = batch)) +
    geom_point(aes(shape = factor(type_labels[group])), size = 5.5, alpha = .7) +
    labs(title = "", x = paste0("PC1: ", pca_var_per[1], "% variance"), 
         y = paste0("PC2: ", pca_var_per[2], "% variance")) +
    scale_shape_manual(values = c(15, 16)) +
    scale_colour_manual(labels = grps_labels, values = colors) +
    theme_bw() +
    theme(strip.text.y = element_text(), 
          axis.text = element_text(colour = "black", size = 15),
          axis.title.x = element_text(colour = "black", face = "bold", size = 15),
          axis.title.y = element_text(colour = "black", face = "bold", size = 15),
          axis.text.x = element_text(colour = "black", size = 15),
          legend.background = element_rect(fill = "white"),
          panel.grid.major = element_line(colour = "white"),
          legend.title = element_text(colour = "black", size = 15),
          legend.text = element_text(colour = "black", size = 15),
          strip.text.x = element_text(colour = "black", size = 15),
          plot.title = element_text(size = 15, face = "bold", hjust = .5, margin = margin(b = 2, unit = "pt")),
          panel.border = element_rect(colour = "black", size = 2),
          legend.key = element_blank()) +
    theme(legend.position = "right")  
  dev.off()
  
  return(pca_var_per)
}

# PCA Before Correction for Training and Validation Sets
run_pca(training$data, training$batch, training$group, type, grps, colors, "before_correction_training")
run_pca(validation$data, validation$batch, validation$group, type, grps, colors, "before_correction_validation")

# Correct batch using ComBat
correct_batch <- function(data, batch, file_prefix) {
  corrected_data <- ComBat(dat = t(data), batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  corrected_data_ <- as.data.frame(cbind(sample = rownames(t(corrected_data)), t(corrected_data)))
  write.table(corrected_data_, paste0(file_prefix, "_corrected.txt"), quote = TRUE, sep = "\t", row.names = TRUE)
  write.csv(corrected_data_, paste0(file_prefix, "_corrected.csv"))
  return(corrected_data_)
}

# Correct batches for training and validation sets
corrected_training <- correct_batch(training$data, training$batch, "corrected_training")
corrected_validation <- correct_batch(validation$data, validation$batch, "corrected_validation")

# PCA After Correction for Training and Validation Sets
run_pca(corrected_training, training$batch, training$group, type, grps, colors, "after_correction_training")
run_pca(corrected_validation, validation$batch, validation$group, type, grps, colors, "after_correction_validation")
