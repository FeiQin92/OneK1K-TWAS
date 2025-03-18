# OneK1K-TWAS
Identification of immune cell type-specific susceptibility genes in multiple cancers using transcriptome-wide association studies

## Author
Fei Qin, Jianxin Shi, Kai Yu

## Description
To identify more susceptibility genes through TWAS, we implemented three distinct elastic net–based approaches to build gene expression prediction models for each cell type: 1) the traditional “targetC” model, which predicts expression using data only from the targeted cell type; 2) the “S+targetC” model, following Thompson et al., which incorporates both shared and cell type-specific expression components for the targeted cell type; and 3) a novel “S+allC” model proposed by us, which integrates shared and cell type-specific expression components across all cell types to further improve prediction performance in the targeted cell type. Prediction accuracy was quantified as R^2, where R is the Pearson correlation coefficient between the predicted and the true expression values, estimated via cross-validations. Ultimately, 6,117~9,058 unique genes with R2 > 0.01 in any of these three models were retained for each of the 14 cell types. These models/genes were used for subsequent TWAS analyses.
