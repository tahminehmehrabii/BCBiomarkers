library(GEOquery); library(dplyr); library(data.table)
library(org.Hs.eg.db); library(annotate); library(preprocessCore)
library(EnsDb.Hsapiens.v79)

setwd("C:/Users/MHR/Desktop/PDAC")

gse_list <- lapply(c("GSE71989", "GSE46234", "GSE15471"), function(gse) getGEO(gse, destdir = getwd())[[1]])

phen_list <- lapply(gse_list, pData)
feat_list <- lapply(gse_list, fData)
expr_list <- lapply(gse_list, function(gse) log2(exprs(gse)))

commonGenes <- Reduce(intersect, lapply(feat_list, function(f) f$`Gene Symbol`))

editSymbol <- function(featData) {
  featData$`Gene Symbol` <- sapply(strsplit(featData$`Gene Symbol`, " /// "), function(x) x[x %in% commonGenes][1])
  featData <- featData[!duplicated(featData$`Gene Symbol`) & featData$`Gene Symbol` != "", c("ID", "Gene Symbol")]
  colnames(featData) <- c("ID", "symbol")
  featData
}

feat_list <- lapply(feat_list, editSymbol)
expr_list <- Map(function(e, f) merge(f, as.data.frame(cbind(ID = rownames(e), e)), by = "ID"), expr_list, feat_list)
expr_list <- lapply(expr_list, function(e) {
  rownames(e) <- e$symbol
  e <- e[, -c(1, 2)]
  na.omit(e)
})

lapply(seq_along(expr_list), function(i) write.csv(expr_list[[i]], paste0("expr", c("15471", "71989", "46234")[i], ".csv"), row.names = TRUE))

defineGroup <- function(phen, patternT, patternN) {
  phen$group <- ifelse(grepl(patternT, phen$source_name_ch1), "T", ifelse(grepl(patternN, phen$source_name_ch1), "N", phen$source_name_ch1))
  write.csv(data.frame(sample = rownames(phen), group = phen$group), paste0("phen", c("15471", "71989", "46234")[which(phen_list == phen)], ".csv"), row.names = FALSE)
}

defineGroup(phen_list[[1]], "T", "N")
defineGroup(phen_list[[2]], "Tissue from PDAC patients", "human normal pancreatic tissue")
defineGroup(phen_list[[3]], "PDAC pancreas", "Normal pancreas")

expr_list <- lapply(seq_along(expr_list), function(i) {
  expr <- as.data.frame(fread(paste0("expr", c("15471", "71989", "46234")[i], ".csv")))
  rownames(expr) <- expr$V1; expr <- expr[, -1]
  expr[Reduce(intersect, lapply(expr_list, rownames)), ]
})
expr_list <- lapply(expr_list, function(e) as.data.frame(t(e)))

phen_list <- lapply(seq_along(phen_list), function(i) {
  phen <- as.data.frame(fread(paste0("phen", c("15471", "71989", "46234")[i], ".csv"), header = TRUE))
  expr <- cbind(sample = rownames(expr_list[[i]]), expr_list[[i]])
  expr <- merge(phen, expr, by = "sample")
  rownames(expr) <- expr$sample; expr[, -1]
})

expr_list <- Map(function(e, b) cbind(batch = b, e), phen_list, 1:3)
mergedExpr <- do.call(rbind, expr_list)
mergedExpr$group <- as.numeric(factor(mergedExpr$group, levels = c("N", "T")))
mergedExpr[, 3:ncol(mergedExpr)] <- lapply(mergedExpr[, 3:ncol(mergedExpr)], as.numeric)
write.csv(mergedExpr, "mergeDatasets.csv")
