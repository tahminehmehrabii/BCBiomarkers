library(data.table)
library(glmnet)
library(ggplot2)
library(dplyr)

setwd("C:/Users/MHR/Desktop/PDAC")

# Load and process train data
trainSet <- fread("corrected_training_set.csv")
trainGroup <- factor(trainSet$group, labels = c("Normal", "Tumor"))
trainData <- trainSet[,-c(1:4)]
keyGenes <- fread("essentialNodes.txt", header = FALSE)$V1
trainData <- trainData[, colnames(trainData) %in% keyGenes, with = FALSE]
train <- cbind(group = trainGroup, trainData)

# LASSO logistic regression
set.seed(123)
x <- model.matrix(group ~ ., train)[,-1]
y <- as.integer(train$group == "Tumor")
cv.lasso <- cv.glmnet(x, y, alpha = 1, family = "binomial", type.measure = "deviance", nfolds = 10)

png("Lambda_BinomialDeviance.png", width = 1800, height = 1800, res = 300)
plot(cv.lasso, cex.lab = 1.4, cex.axis = 1.2)
dev.off()

write.table(cv.lasso$lambda.min, "lambda_min.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Final model and non-zero coefficients
lasso.model <- glmnet(x, y, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)
coef.nonzero <- coef(lasso.model)[-1,1]
df.coef.nonzero <- data.frame(Gene = names(coef.nonzero)[coef.nonzero != 0], Coefficient = coef.nonzero[coef.nonzero != 0])
write.csv(df.coef.nonzero, "coefficients.csv", row.names = FALSE)

write.table(df.coef.nonzero$Gene, "lassoGenes.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Plot non-zero coefficients
df.coef.nonzero <- df.coef.nonzero %>% arrange(Coefficient)
df.coef.nonzero$Gene <- factor(df.coef.nonzero$Gene, levels = df.coef.nonzero$Gene)

png("nonZeroCoefGenes.png", width = 1600, height = 1400, res = 300)
ggplot(df.coef.nonzero, aes(x = Coefficient, y = Gene, fill = Coefficient > 0)) + 
  geom_bar(stat = "identity", width = 0.65) +
  geom_text(aes(label = round(Coefficient, 3)), hjust = ifelse(df.coef.nonzero$Coefficient > 0, -0.1, 1.1), size = 3) +
  theme_bw() +
  xlab("Coefficient") + ylab("") +
  scale_x_continuous(limits = c(min(df.coef.nonzero$Coefficient) - 0.1, max(df.coef.nonzero$Coefficient) + 0.2)) +
  theme(axis.title.x = element_text(size = 13),
        axis.text = element_text(size = 10, colour = "black"),
        legend.position = "none")
dev.off()
