# OneK1K-TWAS
Identification of immune cell type-specific susceptibility genes in multiple cancers using transcriptome-wide association studies

## Author
Fei Qin, Jianxin Shi, Kai Yu

## Description
To identify more susceptibility genes through TWAS, we implemented three distinct elastic net–based approaches to build gene expression prediction models for each cell type: 1) the traditional “targetC” model, which predicts expression using data only from the targeted cell type; 2) the “S+targetC” model, following Thompson et al., which incorporates both shared and cell type-specific expression components for the targeted cell type; and 3) a novel “S+allC” model proposed by us, which integrates shared and cell type-specific expression components across all cell types to further improve prediction performance in the targeted cell type. Prediction accuracy was quantified as R^2, where R is the Pearson correlation coefficient between the predicted and the true expression values, estimated via cross-validations. These models/genes were used for subsequent TWAS analyses.

## Installation
```r
install.packages("devtools")
library(devtools)
install_github("FeiQin92/OneK1K-TWAS")
```

## Generating prediction models
```r
setwd(paste0(".../OneK1K-TWAS-main/data"))
######################################################################################   
############### 1. Fiting hom and het components for each cell type ##################
######################################################################################   

genename <- "ENSG00000000457"
X1 <- paste0("Exp/", genename, "/snpexp")
Y1 <- paste0("Exp/", genename, "/y/")
cov1 <- "Cov/"
Out1 <- "results/"
SNP1 <- paste0("Exp/", genename, "/snppos.txt")
Pre1 <- genename

library(scTWAS)
hom_het_model(X_file=X1,  # Genotype matrix, individuals are row names, no column names. Check examples data for more details.
              Y_file_dir=Y1, # Direction containing only the expression files, individuals are row names, no column names. Check examples data for more details.
              cov_file_dir=cov1, # Cov file direction, should share the same prefix as the expression files Check examples data for more details.
              out_dir=Out1,  # output CV predictors and other coefficients.
              snps=SNP1,  # A tab-delimited txt file containing snp information. Six columns.
              gene_name=Pre1) # Add this prefix to all the saved results files  
```

```r
#########################################################################################
######  2. Fitting models using hom and het components generated in the first step ######
#########################################################################################
Model_fit(Gname=genename)
```

```r
############################################################################################
###### 3. Combing coefficients from above two steps and prepare for weights file   #########
############################################################################################
Weights_pre(Gname=genename, in_dir="results/")
```

```r
#########################################################################################
############# 4. Generating weights files required for FUSION TWAS ######################
#########################################################################################   
Weights_gen(Gname=genename, in_dir="results/", out_dir="weights/")
```
