library(data.table)
setwd("C:/Users/MHR/Desktop/PDAC")

samples <- fread("mergeDatasets.csv")
setDF(samples)
rownames(samples) <- samples[[1]]
samples <- samples[, -1, drop = FALSE]

split_samples <- function(sample_names, ratio = 0.7) {
  train_idx <- sample(seq_along(sample_names), ceiling(ratio * length(sample_names)))
  list(train = sample_names[train_idx], valid = sample_names[-train_idx])
}

cancerSamples <- rownames(samples)[samples$group == 2]
normalSamples <- rownames(samples)[samples$group == 1]

cancerSplit <- split_samples(cancerSamples)
normalSplit <- split_samples(normalSamples)

trainingSamples <- c(cancerSplit$train, normalSplit$train)
validationSamples <- c(cancerSplit$valid, normalSplit$valid)

write.csv(samples[trainingSamples, , drop = FALSE], "training_set.csv")
write.csv(samples[validationSamples, , drop = FALSE], "validation_set.csv")
