# A Bioinformatics and Machine Learning Pipeline for Discovering Diagnostic Gene Biomarkers in Breast Cancer

# Overview
This github repository contains the data files and analysis code used for the scientific Project titled "Integrating machine learning and bioinformatics approaches for identifying novel diagnostic gene biomarkers in Breast cancer". The files are organized into four folders:

Data: which contains all the transcriptomic data required to perform the analyses.
Codes: contains the R scripts to reproduce all analyses.
Results: contains all the results produced by the R scripts.

# Reproducing the results
To reproduce the results, the following R scripts were executed:

1. mergeDatasets.R was run to merge the five microarray datasets (GSE42568, GSE61304) based on their common genes.
2. trainValidSplit.R was executed to split the merged datasets into training and validation datasets.
3. batchCorrection.R was used to correct batch effects in the training and validation datasets.
4. differentialExpressionAnalysis.R was run to identify differentially expressed genes (DEGs) in the training dataset.
5. enrichmentAnalysis.R was performed to gain deeper insight into the biological significance of the DEGs.
6. geneCoexpressionAnalysis.R was executed to construct co-expression modules on the training dataset using the CEMiTool package with default settings.
7. overlappedGenes.R was used to identify overlapping genes between the DEGs and the genes within the most significant module identified by CEMiTool.
8. centralityAnalysis.R was run to identify key genes in the protein–protein interaction (PPI) network.
9. LASSO.R was applied to identify candidate diagnostic genes among the key genes in the training dataset.
10. ROC.R was executed to evaluate the sensitivity and specificity of the candidate diagnostic genes in both training and validation datasets.
11. trainEvaluateML.R was used to verify the ability of the diagnostic genes to distinguish between normal and tumor samples using two machine learning models: Random Forest (RF) and Support Vector Machine (SVM).
12.immuneAnalysis.R was performed to assess immune cell infiltration in breast cancer and to analyze the correlation between biomarkers and immune cells.


# Required software
1. R (4.4.2)
2. RStudio version: 2024.12.0
3. GEOquery (2.74.0)
4. org.Hs.eg.db (3.20.0)
5. EnsDb.Hsapiens.v79 (2.99.0)
6. ggplot2 (3.5.1)
7. sva (3.54.0)
8. igraph (2.1.4)
9. CINNA (1.2.2)
10. ggraph (2.2.1)
11. limma (3.62.2)
12. clusterProfiler (4.14.6)
13. CEMiTool (1.30.0)
14. pheatmap (1.0.12)
15. glmnet (4.1.8)
16. pROC (1.18.5)
17. MLmetrics (1.1.3)
18. randomForest (4.7.1.2)
19. caret (7.0.1)
20. e1071 (1.7.16)
21. data.table (1.17.0)
22. RColorBrewer (1.1.3)
23. preprocessCore (1.68.0)
24. dplyr (1.1.4)
25. annotate (1.84.0)
26. forcats (1.0.0)
27. ggpubr (0.6.0)
28. readr (2.1.5)
29. graphlayouts (1.2.2)
30. enrichplot (1.26.6)
31. AnnotationDbi (1.68.0)
32. ggvenn (0.1.10)
33. ggrepel (0.9.6)
34. neuralnet (1.44.2)
35. nnet (7.3.20)
36. NeuralNetTools (1.5.3)



















