library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(annotate)
library(preprocessCore)
library(EnsDb.Hsapiens.v79)

# Set working directory
setwd("C:/Users/MHR/Desktop/PDAC")

# Load GEO datasets
gse_list <- lapply(c("GSE71989", "GSE46234", "GSE15471"), function(id) getGEO(id, destdir = getwd())[[1]])
names(gse_list) <- c("gse71989", "gse46234", "gse15471")

# Extract phenotype and feature data
phen_list <- lapply(gse_list, pData)
feat_list <- lapply(gse_list, fData)
expr_list <- lapply(gse_list, function(gse) log2(exprs(gse) + 1))  # Added +1 to avoid log(0)

# Find common genes
commonGenes <- Reduce(intersect, lapply(feat_list, function(feat) feat$`Gene Symbol`))

# Function to edit gene symbols
editSymbol <- function(featData) {
  featData$`Gene Symbol` <- sapply(strsplit(featData$`Gene Symbol`, " /// "), function(x) x[x %in% commonGenes][1])
  featData <- featData[!is.na(featData$`Gene Symbol`), c("ID", "Gene Symbol")]
  featData <- featData[!duplicated(featData$`Gene Symbol`), ]  # Remove duplicates
  return(featData)
}

# Apply function to feature data
feat_list <- lapply(feat_list, editSymbol)
names(feat_list) <- names(gse_list)

# Merge expression data with feature data
expr_list <- Map(function(feat, expr) {
  expr <- cbind(ID = rownames(expr), expr) %>% as.data.frame()
  expr <- merge(feat, expr, by = "ID")
  rownames(expr) <- expr$`Gene Symbol`
  expr <- expr[, -c(1,2)]  # Remove ID and Gene Symbol columns
  return(expr)
}, feat_list, expr_list)

# Save expression data
lapply(names(expr_list), function(name) write.csv(expr_list[[name]], paste0(name, ".csv"), row.names = TRUE))

# Identify biological classes
phen_list[["gse15471"]]$group <- ifelse(startsWith(phen_list[["gse15471"]]$title, "T"), "T", "N")
phen_list[["gse46234"]]$group <- ifelse(grepl("Normal pancreas", phen_list[["gse46234"]]$source_name_ch1), "N", "T")

# Save phenotype data
lapply(names(phen_list), function(name) write.csv(data.frame(sample = rownames(phen_list[[name]]), group = phen_list[[name]]$group), paste0(name, "_phen.csv"), row.names = FALSE))

# Load expression and phenotype data
expr_list <- lapply(names(expr_list), function(name) {
  df <- fread(paste0(name, ".csv"))
  rownames(df) <- df$V1
  df <- df[, -1]
  return(df)
})

phen_list <- lapply(names(phen_list), function(name) fread(paste0(name, "_phen.csv")))

# Find intersecting genes
intersectSymbol <- Reduce(intersect, lapply(expr_list, rownames))
expr_list <- lapply(expr_list, function(df) df[intersectSymbol, ])
expr_list <- lapply(expr_list, function(df) as.data.frame(t(df)))

# Merge datasets
expr_list <- Map(function(expr, phen, batch) {
  expr <- cbind(sample = rownames(expr), expr)
  expr <- merge(phen, expr, by = "sample")
  expr$batch <- batch
  return(expr)
}, expr_list, phen_list, 1:3)

mergedExpr <- do.call(rbind, expr_list)
mergedExpr$group <- as.numeric(as.factor(mergedExpr$group))
mergedExpr[, -c(1,2)] <- lapply(mergedExpr[, -c(1,2)], as.numeric)

write.csv(mergedExpr, "mergeDatasets.csv", row.names = FALSE)
