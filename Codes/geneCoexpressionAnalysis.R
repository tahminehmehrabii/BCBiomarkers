library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("CEMiTool")
library(CEMiTool)

# Set working directory
setwd("C:/Users/MHR/Desktop/PDAC")

#### Load and process data ####
trainingSet <- fread("corrected_training_set.csv")
exprData <- t(trainingSet[, -c(1:4)])
phenData <- gsub(1, "Normal", gsub(2, "Tumor", trainingSet$group))
rownames(exprData) <- trainingSet$sample
write.csv(exprData, paste0("GeneExpressionProfileOf", ncol(exprData), "Samples.csv"), quote = F)
write.csv(data.frame(ID = colnames(exprData), Group = phenData), 
          paste0("DiseasePhenotypesOf", ncol(exprData), "Samples.csv"), quote = F)

#### WGCNA using CEMiTool ####
cem <- cemitool(expr = exprData, annot = data.frame(ID = colnames(exprData), Group = phenData),
                gmt = read_gmt(system.file("extdata", "pathways.gmt", package = "CEMiTool")),
                interactions = read.delim(system.file("extdata", "interactions.tsv", package = "CEMiTool")),
                sample_name_column = "ID", class_column = "Group", plot = TRUE, verbose = TRUE)

# Generate report, tables, and plots
generate_report(cem, directory = "geneCoexpressionAnalysis/Report", force = TRUE, output_format = "html_document")
write_files(cem, directory = "geneCoexpressionAnalysis/Tables", force = TRUE)
save_plots(cem, value = c("all"), force = TRUE, directory = "geneCoexpressionAnalysis/Plots")

# Save specific plots
plots <- c("beta_r2", "mean_k", "gsea", "ora_M1")
for (plot in plots) {
  png(paste0("geneCoexpressionAnalysis/", plot, ".png"), height = 1600, width = 1600, res = 300)
  show_plot(cem, plot)
  dev.off()
}

# Save significant genes
sigGenes <- fread("geneCoexpressionAnalysis/Tables/module.tsv")$genes
write.table(sigGenes[which(fread("geneCoexpressionAnalysis/Tables/module.tsv")$modules == "M1")], 
            "geneCoexpressionAnalysis/modulesGenes.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
