# PDAC-Biomarkers
Overview
This github repository contains the data files and analysis code used for the scientific paper titled "Integrating machine learning and bioinformatics approaches for identifying novel diagnostic gene biomarkers in Breast cancer". 

Data: which contains all the transcriptomic data required to perform the analyses.
Codes: contains the R scripts to reproduce all analyses.
Results: contains all the results produced by the R scripts.
Reproducing the results
This repository contains all the code necessary to reproduce the results.

git clone https://github.com/ayoub-vaziri/CRCBiomarkers.git path/to/directory

mergeDatasets.R was executed to merge the microarray datasets based on their common genes.

trainValidSplit.R was executed to split the merged datasets into training and validation datasets.

batchCorrection.R was executed to correct batch effects in the training and validation datasets.

differentialExpressionAnalysis.R was executed to identify differentially expressed genes (DEGs) in the training dataset.

enrichmentAnalysis.R was executed to gain deeper insights into the biological significance of the DEGs.

geneCoexpressionAnalysis.R was executed to construct co-expression modules in the training dataset using the CEMiTool package with default settings.

overlappedGenes.R was executed to identify genes overlapping between the DEGs and the most significant module identified by CEMiTool.

centralityAnalysis.R was executed to identify key genes in the protein-protein interaction (PPI) network.

LASSO.R was executed to identify candidate diagnostic genes from the set of key genes in the training dataset.

ROC.R was executed to evaluate the sensitivity and specificity of candidate diagnostic genes in both the training and validation datasets.

trainEvaluateML.R was executed to verify the accuracy of diagnostic genes in distinguishing between normal and tumor samples using machine learning models.

SignallingAnalysis.R was executed to assess signaling pathways in breast cancer (BC) and the correlation between biomarkers and signaling pathways.
