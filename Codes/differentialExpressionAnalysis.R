library(limma)
library(readr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# Set working directory
setwd("C:/Users/MHR/Desktop/Breast")

#### function log2 transform ####
log2trans <- function(expr) {
  quan <- as.numeric(quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (quan[5] > 100) || (quan[6] - quan[1] > 50 && quan[2] > 0)
  if (LogC) {
    expr[which(expr <= 0)] <- NaN
    expr <- log2(expr)
  }
  return(expr)
}

#### DEGs function ####
degs <- function(expr, group) {
  design.mat <- model.matrix(~ group + 0)
  colnames(design.mat) <- levels(group)
  lmfit <- lmFit(expr, design.mat)
  cont.mat <- makeContrasts(contrasts = "Tumor-Normal", levels = design.mat)
  cont.fit <- contrasts.fit(lmfit, cont.mat)
  fit <- eBayes(cont.fit, 0.01)
  toptable <- topTable(fit = fit, number = Inf, adjust.method = "fdr", sort.by = "B")
  return(toptable)
}

#### Reading data and preparing ####
trainingSet <- as.data.frame(fread("corrected_training_set.csv"))

group <- trainingSet$group
group <- gsub(1, "N", group)
group <- gsub(2, "C", group)
group <- as.factor(group)
levels(group) <- c("Tumor", "Normal")

rownames(trainingSet) <- trainingSet$V1
trainingSet <- trainingSet[, -c(1, 2)]
trainingSet <- as.data.frame(t(trainingSet))

trainingSet <- log2trans(trainingSet)
tT <- degs(trainingSet, group)

tT$Gene.symbol <- rownames(tT)
tT <- tT[, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]
rownames(tT) <- NULL

write_tsv(tT, "topTable.tsv")
write_csv(tT, "topTable.csv")

#### Labeling step for Volcano Plot ####
tT$threshold <- "Not significant"
tT$threshold[tT$logFC > 1 & tT$adj.P.Val < 0.05] <- "Up-regulated"
tT$threshold[tT$logFC < -1 & tT$adj.P.Val < 0.05] <- "Down-regulated"

# Select Top 10 Up and Down genes based on lowest adj.P.Val
top_up <- tT[tT$threshold == "Up-regulated", ]
top_up <- top_up[order(top_up$adj.P.Val), ][1:min(10, nrow(top_up)), ]

top_down <- tT[tT$threshold == "Down-regulated", ]
top_down <- top_down[order(top_down$adj.P.Val), ][1:min(10, nrow(top_down)), ]

top_genes <- rbind(top_up, top_down)

#### رسم Volcano Plot ####
# Create a volcano plot with a white background and black labels
volcano_plot <- ggplot(tT, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = "blue", 
                                "Not significant" = "gray", 
                                "Up-regulated" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot (Top 10 Genes Labeled for Each Direction)",
       x = "log2 Fold Change",
       y = "-log10(P-value)") +
  theme_bw() +  
  
# Tagging up-regulated genes
  geom_text_repel(data = tT[tT$logFC > 1 & tT$adj.P.Val < 0.05, ][1:10, ],
                  aes(label = Gene.symbol),
                  size = 3,
                  max.overlaps = 100,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  color = "black") +  # رنگ برچسب‌ها مشکی
  
# Tagging down-regulated genes
  geom_text_repel(data = tT[tT$logFC < -1 & tT$adj.P.Val < 0.05, ][1:10, ],
                  aes(label = Gene.symbol),
                  size = 3,
                  max.overlaps = 100,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.color = "grey50",
                  color = "black")    # رنگ برچسب‌ها مشکی

# Show plot
print(volcano_plot)

# Save as high quality PNG
ggsave("volcano_plot_labeled.png",
       plot = volcano_plot,
       width = 7,
       height = 6,
       dpi = 300)

# Assume that the data matrix and metadata are ready
# Data: trainingSet (expression matrix with genes in rows and samples in columns)
# Group information: group_info (data frame with group information for each sample)

# Define groups for samples

group_info <- data.frame(Group = group)
rownames(group_info) <- colnames(trainingSet)

# Define colors for groups
ann_colors <- list(
  Group = c(Normal = "green", Tumor = "red")
)

# Draw a heatmap and save it as a PNG file
heatmap_plot <- pheatmap(trainingSet,
                         annotation_col = group_info,
                         annotation_row = group_info,  # در صورت تمایل
                         annotation_colors = ann_colors,
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_rownames = FALSE,
                         show_colnames = FALSE,
                         fontsize = 10,
                         border_color = NA)

# Save heatmap as high quality PNG file
ggsave("heatmap_plot.png",
       plot = heatmap_plot,
       width = 8,
       height = 6,
       dpi = 300)

#########

##########
# Loading libraries
library(data.table)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Upload file
tT <- as.data.frame(fread("C:/Users/MHR/Desktop/Breast/corrected_training_set.csv"))

# If Gene.symbol exists, set it as rownames
rownames(tT) <- tT$Gene.symbol

# Filtering differentially expressed genes (DEGs)
DEGs <- tT %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC)))

# Extract expression data only from columns related to examples (let's say from column 4 onwards)
expr_data <- tT[DEGs$Gene.symbol, 4:ncol(tT)]

# Define annotation for group of samples (assuming group variable is already defined)
# For example, group <- c("Tumor", "Tumor", "Normal", ...) as many columns as expr_data

annotation_col <- data.frame(Group = factor(group))
rownames(annotation_col) <- colnames(expr_data)

# Define annotation for rows (genes)
annotation_row <- data.frame(Regulated = ifelse(DEGs$logFC > 1, "up", "down"))
rownames(annotation_row) <- rownames(expr_data)

# Colors for annotation
ann_colors <- list(
  Group = c(Normal = "#00BCD4", Tumor = "#E1BEE7"),
  Regulated = c(up = "red", down = "green")
)

# Color palette for heatmap
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Draw a heatmap
pheatmap(expr_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = heat_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize = 10,
         border_color = NA,
         main = paste("Heatmap of", nrow(expr_data), "DEGs"))

# Save heatmap
png("Heatmap_DEGs.png", width = 1200, height = 1000, res = 150)
pheatmap(expr_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = heat_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize = 10,
         border_color = NA,
         main = paste("Heatmap of", nrow(expr_data), "DEGs"))
dev.off()



################
### Loading libraries
library(limma)
library(readr)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)

### Specify the file path
setwd("C:/Users/MHR/Desktop/Breast")

### log2 transform function
log2trans <- function(expr) {
  quan <- as.numeric(quantile(expr, c(0.0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  LogC <- (quan[5] > 100) || (quan[6] - quan[1] > 50 && quan[2] > 0)
  if (LogC) {
    expr[which(expr <= 0)] <- NaN
    expr <- log2(expr)
  }
  return(expr)
}

### DEGs function
degs <- function(expr, group) {
  design.mat <- model.matrix(~ group + 0)
  colnames(design.mat) <- levels(group)
  lmfit <- lmFit(expr, design.mat)
  cont.mat <- makeContrasts(contrasts = "Tumor-Normal", levels = design.mat)
  cont.fit <- contrasts.fit(lmfit, cont.mat)
  fit <- eBayes(cont.fit, 0.01)
  toptable <- topTable(fit = fit, number = Inf, adjust.method = "fdr", sort.by = "B")
  return(toptable)
}

### Reading and preparing data
raw_data <- as.data.frame(fread("corrected_training_set.csv"))
group <- raw_data$group
group <- gsub(1, "N", group)
group <- gsub(2, "C", group)
group <- as.factor(group)
levels(group) <- c("Tumor", "Normal")

rownames(raw_data) <- raw_data$V1
expr_data <- raw_data[, -c(1, 2)]  
expr_data <- as.data.frame(t(expr_data))  
expr_data <- log2trans(expr_data)

### Identification of DEGs genes
tT <- degs(expr_data, group)
tT$Gene.symbol <- rownames(tT)
tT <- tT[, c("Gene.symbol", "logFC", "P.Value", "adj.P.Val")]

### Save output DEGs
write_csv(tT, "topTable.csv")

### Extracting differentially expressed genes (DEGs) with logFC > 1 and adj.P.Val < 0.05
DEGs <- tT %>%
  filter(abs(logFC) > 1 & adj.P.Val < 0.05) %>%
  arrange(desc(abs(logFC)))

### Selecting data related to DEGs genes
expr_DEGs <- expr_data[DEGs$Gene.symbol, ]

### Sorting samples by group (Normal and Tumor)
sample_order <- order(group)  # اول Normal، بعد Tumor
expr_DEGs <- expr_DEGs[, sample_order]

### Annotation definition for examples
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr_DEGs)
annotation_col <- annotation_col[sample_order, , drop = FALSE]

### Definition of annotation for genes
annotation_row <- data.frame(Regulated = ifelse(DEGs$logFC > 1, "up", "down"))
rownames(annotation_row) <- rownames(expr_DEGs)

### Definition of colors
ann_colors <- list(
  Group = c(Normal = "#00BCD4", Tumor = "#E1BEE7"),
  Regulated = c(up = "red", down = "green")
)

### Draw and save the Heatmap
heat_colors <- colorRampPalette(c("blue", "white", "red"))(100)
png("Heatmap_DEGs.png", width = 1200, height = 1000, res = 150)
pheatmap(expr_DEGs,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = ann_colors,
         color = heat_colors,
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize = 10,
         border_color = NA,
         main = paste("Heatmap of", nrow(expr_DEGs), "DEGs"))
dev.off()

### Making Volcano Plot
tT$threshold <- "Not significant"
tT$threshold[tT$logFC > 1 & tT$adj.P.Val < 0.05] <- "Up-regulated"
tT$threshold[tT$logFC < -1 & tT$adj.P.Val < 0.05] <- "Down-regulated"

top_up <- tT[tT$threshold == "Up-regulated", ] %>% arrange(adj.P.Val) %>% head(10)
top_down <- tT[tT$threshold == "Down-regulated", ] %>% arrange(adj.P.Val) %>% head(10)

volcano_plot <- ggplot(tT, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Down-regulated" = "blue", "Not significant" = "gray", "Up-regulated" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(title = "Volcano Plot (Top 10 Genes Labeled for Each Direction)",
       x = "log2 Fold Change", y = "-log10(P-value)") +
  theme_bw() +
  geom_text_repel(data = top_up, aes(label = Gene.symbol), size = 3, color = "black") +
  geom_text_repel(data = top_down, aes(label = Gene.symbol), size = 3, color = "black")

ggsave("volcano_plot_labeled.png", plot = volcano_plot, width = 7, height = 6, dpi = 300)

# Sort samples by group (Normal and Tumor)
sample_order <- order(group) 
expr_DEGs <- expr_DEGs[, sample_order]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

#############
