# PDAC-Biomarkers
Overview
This github repository contains the data files and analysis code used for the scientific paper titled "Integrating machine learning and bioinformatics approaches for identifying novel diagnostic gene biomarkers in Breast cancer". 

Data: which contains all the transcriptomic data required to perform the analyses.
Codes: contains the R scripts to reproduce all analyses.
Results: contains all the results produced by the R scripts.
Reproducing the results
This repository contains all the code necessary to reproduce the results.

git clone (https://github.com/tahminehmehrabii/Breast Cancer-Biomarkers.git)

1. mergeDatasets.R was executed to merge the microarray datasets based on their common genes.

2. trainValidSplit.R was executed to split the merged datasets into training and validation datasets.

3. batchCorrection.R was executed to correct batch effects in the training and validation datasets.

4. differentialExpressionAnalysis.R was executed to identify differentially expressed genes (DEGs) in the training dataset.

5. enrichmentAnalysis.R was executed to gain deeper insights into the biological significance of the DEGs.

6. geneCoexpressionAnalysis.R was executed to construct co-expression modules in the training dataset using the CEMiTool package with default settings.

7. overlappedGenes.R was executed to identify genes overlapping between the DEGs and the most significant module identified by CEMiTool.

8. centralityAnalysis.R was executed to identify key genes in the protein-protein interaction (PPI) network.

9. LASSO.R was executed to identify candidate diagnostic genes from the set of key genes in the training dataset.

10. ROC.R was executed to evaluate the sensitivity and specificity of candidate diagnostic genes in both the training and validation datasets.

11. trainEvaluateML.R was executed to verify the accuracy of diagnostic genes in distinguishing between normal and tumor samples using machine learning models.

12. SignallingAnalysis.R was executed to assess signaling pathways in breast cancer (BC) and the correlation between biomarkers and signaling pathways.
