######################
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(annotate)
library(preprocessCore)
library(EnsDb.Hsapiens.v79)

# Set working directory
setwd("C:/Users/MHR/Desktop/PDAC")

#### Load data ####
gse71989 <- getGEO("GSE71989", destdir = getwd())[[1]]
gse46234 <- getGEO("GSE46234", destdir = getwd())[[1]]
gse15471 <- getGEO("GSE15471", destdir = getwd())[[1]]

# Extract phenotype and feature data
phen71989 <- pData(gse71989)
phen46234 <- pData(gse46234)
phen15471 <- pData(gse15471)

feat71989 <- fData(gse71989)
feat46234 <- fData(gse46234)
feat15471 <- fData(gse15471)

# Expression data
expr15471 <- exprs(gse15471)
expr71989 <- log2(exprs(gse71989))
expr46234 <- log2(exprs(gse46234))



# Find common genes
commonGenes <- Reduce(intersect, list(feat71989$`Gene Symbol`, 
                                      feat46234$`Gene Symbol`, 
                                      feat15471$`Gene Symbol`))

# Function to edit symbol names
editSymbol <- function(featData, geneColName) {
  rows <- which(unlist(lapply(strsplit(featData[,geneColName], " /// "), function(x) length(x) > 1)))
  for (x in rows) {
    ss <- strsplit(featData[,geneColName][x], " /// ")[[1]]
    idx <- which(ss %in% commonGenes)
    if (length(idx) > 0) {
      featData[,geneColName][x] <- ss[idx[1]]
    }
  }
  return(featData)
}

# Apply function and select correct columns
feat15471 <- editSymbol(feat15471, "Gene Symbol")
feat46234 <- editSymbol(feat46234, "Gene Symbol")
feat71989 <- editSymbol(feat71989, "Gene Symbol")

feat15471 <- feat15471[,c("ID", "Gene Symbol")]
feat46234 <- feat46234[,c("ID", "Gene Symbol")]
feat71989 <- feat71989[,c("ID", "Gene Symbol")]

colnames(feat15471) <- c("ID", "symbol")
colnames(feat46234) <- c("ID", "symbol")
colnames(feat71989) <- c("ID", "symbol")

# Convert expression data to data.frame and merge with feature data
expr15471 <- cbind(ID = rownames(expr15471), expr15471) %>% as.data.frame()
expr46234 <- cbind(ID = rownames(expr46234), expr46234) %>% as.data.frame()
expr71989 <- cbind(ID = rownames(expr71989), expr71989) %>% as.data.frame()

expr15471 <- merge(feat15471, expr15471, by = "ID")
expr46234 <- merge(feat46234, expr46234, by = "ID")
expr71989 <- merge(feat71989, expr71989, by = "ID")

# Remove empty symbols
expr15471 <- expr15471[which(expr15471$symbol != ""),]
expr46234 <- expr46234[which(expr46234$symbol != ""),]
expr71989 <- expr71989[which(expr71989$symbol != ""),]

# Remove duplicate symbols
expr15471 <- expr15471[!duplicated(expr15471$symbol),]
expr46234 <- expr46234[!duplicated(expr46234$symbol),]
expr71989 <- expr71989[!duplicated(expr71989$symbol),]

# Set row names and remove redundant columns
rownames(expr15471) <- expr15471$symbol
rownames(expr46234) <- expr46234$symbol
rownames(expr71989) <- expr71989$symbol

expr15471 <- expr15471[,-c(1,2)]
expr46234 <- expr46234[,-c(1,2)]
expr71989 <- expr71989[,-c(1,2)]

# Remove NA values
expr15471 <- na.omit(expr15471)
expr46234 <- na.omit(expr46234)
expr71989 <- na.omit(expr71989)


# Save raw expression data
write.csv(expr15471, "expr15471.csv", row.names = TRUE)
write.csv(expr71989, "expr71989.csv", row.names = TRUE)
write.csv(expr46234, "expr46234.csv", row.names = TRUE)


#### Identifying biological classes ####
########################################
phen15471$title <- ifelse(startsWith(phen15471$title, "T"), "T", "N")
phen15471_ <- data.frame(sample = rownames(phen15471), group = phen15471$title)
write.csv(phen15471_, "phen15471.csv", row.names = FALSE)


# ؟؟؟؟؟
phen71989 <- read.csv("phen71989.csv")


phen46234$source_name_ch1 <- ifelse(grepl("Normal pancreas", phen46234$source_name_ch1), "N",
                                    ifelse(grepl("PDAC pancreas", phen46234$source_name_ch1), "T", phen46234$source_name_ch1))
phen46234_ <- data.frame(sample = rownames(phen46234), group = phen46234$source_name_ch1)
write.csv(phen46234_, "phen46234.csv", row.names = FALSE)



#### Load expression data ####
##############################
expr15471 <- as.data.frame(fread("expr15471.csv"))
rownames(expr15471) <- expr15471$V1
expr15471 <- expr15471[,-1]

expr46234 <- as.data.frame(fread("expr46234.csv"))
rownames(expr46234) <- expr46234$V1
expr46234 <- expr46234[,-1]

expr71989 <- as.data.frame(fread("expr71989.csv"))
rownames(expr71989) <- expr71989$V1
expr71989 <- expr71989[,-1]


intersectSymbol <- Reduce(intersect, 
                          list(
                            rownames(expr15471),
                            rownames(expr46234),
                            rownames(expr71989)))

expr15471 <- expr15471[intersectSymbol,]
expr46234 <- expr46234[intersectSymbol,]
expr71989 <- expr71989[intersectSymbol,]

expr15471 <- t(expr15471) %>% as.data.frame()
expr46234 <- t(expr46234) %>% as.data.frame()
expr71989 <- t(expr71989) %>% as.data.frame()



#### Load phenotypic data ####
##############################
phen15471 <- as.data.frame(fread("phen15471.csv", header = TRUE))
expr15471 <- cbind(sample = rownames(expr15471), expr15471)
expr15471 <- merge(phen15471, expr15471, by="sample")
rownames(expr15471) <- expr15471$sample
expr15471 <- expr15471[,-1]

phen46234 <- as.data.frame(fread("phen46234.csv", header = TRUE))
expr46234 <- cbind(sample = rownames(expr46234), expr46234)
expr46234 <- merge(phen46234, expr46234, by="sample")
rownames(expr46234) <- expr46234$sample
expr46234 <- expr46234[,-1]

phen71989 <- as.data.frame(fread("phen71989.csv", header = TRUE))
expr71989 <- cbind(sample = rownames(expr71989), expr71989)
expr71989 <- merge(phen71989, expr71989, by = "sample")
rownames(expr71989) <- expr71989$sample
expr71989 <- expr71989[, -1]

#### Merging datasets ####
##########################
expr15471 <- cbind(batch=rep(1, nrow(expr15471)), expr15471)
expr46234<- cbind(batch=rep(2, nrow(expr46234)), expr46234)
expr71989 <- cbind(batch=rep(3, nrow(expr71989)), expr71989)

mergedExpr <- rbind(expr15471, expr46234, expr71989)

mergedExpr$group <- as.factor(mergedExpr$group)
levels(mergedExpr$group) <- c(1,2)
mergedExpr$group <- as.numeric(mergedExpr$group)

for(x in 3:ncol(mergedExpr)) mergedExpr[,x] <- as.numeric(mergedExpr[,x])

write.csv(mergedExpr, "mergeDatasets.csv")
