library(data.table)

setwd("C:/Users/MHR/Desktop/PDAC")

samples <- fread("mergeDatasets.csv")
samples <- samples[, -1, with = FALSE]
rownames(samples) <- samples$V1

# Sample assignment
split_samples <- function(group) {
  samples_group <- rownames(samples)[samples$group == group]
  training <- sample(samples_group, ceiling(0.7 * length(samples_group)))
  validation <- setdiff(samples_group, training)
  return(list(training = training, validation = validation))
}

cancerSamples <- split_samples(2)
normalSamples <- split_samples(1)

# Combine training and validation sets
trainingSamples <- union(cancerSamples$training, normalSamples$training)
validationSamples <- union(cancerSamples$validation, normalSamples$validation)

# Create training and validation sets
trainingSet <- samples[rownames(samples) %in% trainingSamples, ]
validationSet <- samples[rownames(samples) %in% validationSamples, ]

# Save results
write.csv(trainingSet, "training_set.csv")
write.csv(validationSet, "validation_set.csv")
