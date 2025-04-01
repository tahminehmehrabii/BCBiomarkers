
######################
# Micro.Array Breast cancer 
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(annotate)
library(preprocessCore)
library(EnsDb.Hsapiens.v79)

# Set the current working directory to the project path
setwd("C:/Users/MHR/Desktop/Breast")

#### Load data ####
###################
gse42568 <- getGEO("GSE42568", destdir = "C:/Users/MHR/Desktop/Breast")[[1]]
gse61304 <- getGEO("GSE61304", destdir = "C:/Users/MHR/Desktop/Breast")[[1]]

phen42568 <- pData(gse42568)
phen61304 <- pData(gse61304)

feat42568 <- fData(gse42568)
feat61304<- fData(gse61304)

# Expression data
expr42568 <- exprs(gse42568)
expr61304 <- log2(exprs(gse61304))
expr61304 <- na.omit(expr61304)


# Find common genes
commonGenes <- Reduce(intersect, list(feat42568$`Gene Symbol`, 
                                      feat61304$`Gene Symbol`))

editSymbol <- function(featData, geneColName) {
  rows <- which(unlist(lapply(strsplit(featData[,geneColName], " /// "), function(x) ifelse(length(x) > 1, TRUE, FALSE))))
  for(x in rows) {
    ss <- strsplit(featData[,geneColName][x], " /// ")[[1]]
    idx <- which(ss %in% commonGenes)[1]
    featData[,geneColName][x] <- ss[idx[1]]
  }
  return(featData)
}

colnames(feat42568)
colnames(feat61304)


feat42568 <- feat42568[,c(1,11)]
colnames(feat42568) <- c("ID", "symbol")

feat61304 <- feat61304[,c(1,11)]
colnames(feat61304) <- c("ID", "symbol")


expr42568 <- cbind(ID=rownames(expr42568), expr42568) %>% as.data.frame()
expr61304 <- cbind(ID=rownames(expr61304), expr61304) %>% as.data.frame()


expr42568  <- merge(feat42568 , expr42568 , by = "ID")
expr61304 <- merge(feat61304, expr61304, by = "ID")

expr42568 <- expr42568[which(expr42568$symbol != ""),]
expr61304 <- expr61304[which(expr61304$symbol != ""),]


expr42568  <- expr42568 [!duplicated(expr42568 $symbol),]
expr61304 <- expr61304[!duplicated(expr61304$symbol),]

rownames(expr42568) <- expr42568$symbol
rownames(expr61304) <- expr61304$symbol


expr42568  <- expr42568 [,-c(1,2)]
expr61304 <- expr61304[, -c(1,2)]


expr42568 <- na.omit(expr42568)
expr61304  <- na.omit(expr61304)

write.csv(expr42568, "expr42568.csv")
write.csv(expr61304, "expr61304.csv")

#### Identifying biological classes ####
########################################

for(x in 1:nrow(phen42568)) {
  if(grepl("normal", phen42568$source_name_ch1[x], ignore.case = TRUE)) {
    phen42568$source_name_ch1[x] <- "N"
  } else if(grepl("cancer", phen42568$source_name_ch1[x], ignore.case = TRUE)) {
    phen42568$source_name_ch1[x] <- "T"
  }
}


for(x in 1:nrow(phen61304)) {
  if(grepl("normal", phen61304$source_name_ch1[x], ignore.case = TRUE)) {
    phen61304$source_name_ch1[x] <- "N"  
  } else if(grepl("tumor", phen61304$source_name_ch1[x], ignore.case = TRUE)) {
    phen61304$source_name_ch1[x] <- "T"  
  }
}


phen42568_ <- as.data.frame(cbind(sample=rownames(phen42568), group=phen42568[,8]))
write.csv(phen42568_, "phen42568.csv")

phen61304_ <- as.data.frame(cbind(sample=rownames(phen61304), group=phen61304[,8]))
write.csv(phen61304_, "phen61304.csv")


#### Load expression data ####
##############################

expr42568 <- as.data.frame(fread("expr42568.csv"))
rownames(expr42568) <- expr42568$V1
expr42568 <- expr42568[,-1]

expr61304 <- as.data.frame(fread("expr61304.csv"))
rownames(expr61304) <- expr61304$V1
expr61304 <- expr61304[,-1]



intersectSymbol <- Reduce(intersect, 
                          list(
                            rownames(expr42568),
                            rownames(expr61304)))

expr42568<- expr42568[intersectSymbol,]
expr61304 <- expr61304[intersectSymbol,]


expr42568 <- t(expr42568) %>% as.data.frame()
expr61304 <- t(expr61304) %>% as.data.frame()


#### Load phenotypic data ####
##############################
phen42568<- as.data.frame(fread("phen42568.csv", header = TRUE))
phen42568 <- phen42568[,-1]
expr42568 <- cbind(sample=rownames(expr42568), expr42568)
expr42568 <- merge(phen42568, expr42568, by="sample")
rownames(expr42568) <- expr42568$sample
expr42568 <- expr42568[,-1]

phen61304<- as.data.frame(fread("phen61304.csv", header = TRUE))
phen61304 <- phen61304[,-1]
expr61304<- cbind(sample=rownames(expr61304), expr61304)
expr61304 <- merge(phen61304, expr61304, by="sample")
rownames(expr61304) <- expr61304$sample
expr61304 <- expr61304[,-1]


#### Merging datasets ####
##########################
expr42568<- cbind(batch=rep(1, nrow(expr42568)), expr42568)
expr61304 <- cbind(batch=rep(2, nrow(expr61304)), expr61304)


mergedExpr <- rbind(expr42568, expr61304)
mergedExpr$group <- as.factor(mergedExpr$group)
levels(mergedExpr$group) <- c(1,2)
mergedExpr$group <- as.numeric(mergedExpr$group)

for(x in 3:ncol(mergedExpr)) mergedExpr[,x] <- as.numeric(mergedExpr[,x])

write.csv(mergedExpr, "mergeDatasets.csv")

