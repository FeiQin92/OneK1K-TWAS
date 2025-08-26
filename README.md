# OneK1K-TWAS
Identification of immune cell type-specific susceptibility genes in multiple cancers using transcriptome-wide association studies

## Author
Fei Qin, Jianxin Shi, Kai Yu

## 1. Introduction to our TWAS framework
We identified susceptibility genes for various cancers using TWAS based on a scRNA-seq dataset comprising 14 immune cell types, along with GWAS summary statistics data for different cancer types. To enhance the power of TWAS, we developed a model to generate genetic predictors (cis-SNPs located ±1 Mb of the gene) for gene expression by incorporating both shared and specific effects across cell types. First, pseudo-bulk expression data for each subject were calculated by aggregating cells within each cell type. Expression values for each gene were then adjusted and residualized over relevant covariates before being decomposed into shared and cell type-specific components. Thereafter, a two-step procedure was implemented to generate an expression prediction model for each gene in each cell type. In the first step, separate models were built using cis-SNPs to predict the shared and cell type-specific gene expression components across all 14 immune cell types. This was achieved through the elastic net, a regularized regression method that balances model complexity and predictive accuracy. In the second step, gene expression within each targeted cell type was predicted using models integrating the predicted shared and cell type-specific expression components from the first step. Finally, weight vectors for each gene and cell type were derived from these models and applied in TWAS to evaluate associations between predicted gene expression and cancer risk across various cell types.

## 2. Installation
```r
#install.packages("devtools")
library(devtools)
install_github("FeiQin92/OneK1K-TWAS")
```

## 3. Download example data for building prediction models
Download the example data from
*[Example_data.zip](https://www.dropbox.com/scl/fo/iksauut206yrep2ouaz9a/APvZfl_XNHT665xJTlHPJw4?rlkey=hz2kc8p93h2jlhk2nplp32h6t&st=5nuvkq2i&dl=1)*. 
```
--Cov
   --genexp_*: Covariates to adjust in each cell type, included sex, age, the top six principal component (PC) scores derived from host genotypic data, and the top two probabilistic estimation of expression residuals (PEER) factors calculated from the expression data.
--Exp
   --ENSG00000000457
       --snpexp: Genotype matrix. Each row is a individual, each column is a cis-SNP for gene ENSG00000000457 (located ±1 Mb of the gene).
       --snppos: A tab-delimited txt file containing snp information with six columns: 1) Chromosome; 2) SNP ID (e.g., 1_168818869); 3) location CM (not exactly necessary for TWAS) 4) Location on chromosome; 5) Allele 1; 6) Allele 2.
       --y
           --genexp_*: Pseudo-bulk expression value for gene ENSG00000000457 in each cell type.
--results: Folder storing temporary files.
--weights: Weights files utlized for TWAS analysis.
```

## 4. Generating prediction models
More details about how to generate prediction models and conducting TWAS have been provided in 
*[vignettes]{}*.


## 5. Predictions models available for all 14 cell types

To identify more susceptibility genes through TWAS, we implemented three distinct elastic net-based approaches to build gene expression prediction models for each cell type: 1) the traditional 'targetC' model, which predicts expression using data only from the targeted cell type; 2) the 'S+targetC' model, following Thompson et al., which incorporates both shared and cell type-specific expression components for the targeted cell type; and 3) a novel 'S+allC' model proposed by us, which integrates shared and cell type-specific expression components across all cell types to further improve prediction performance in the targeted cell type. *[Prediction models](https://www.dropbox.com/scl/fi/7bou9lku437k47uqcma5t/weights_org.tar?rlkey=4ta9pxzzgxs7te1bqx9k7jbvm&st=ry2v6tnb&dl=1)* for all 14 cell types have been provided according to the the format requirements of *[FUSION](https://github.com/gusevlab/fusion_twas)*.

14 Cell type: "CD4_ET", "NK", "CD4_NC", "CD8_S100B", "CD8_ET", "B_IN", "CD8_NC", "B_Mem", "NK_R", "Mono_NC", "Mono_C", "DC", "Plasma", "CD4_SOX4".\

Models were compared (see table below) based on the number of genes selected in each model if prediction performance R2 > 0.01. Finally, genes with prediction performance R2 > 0.01 in either of these three models were retained for further TWAS analyses (“Pooled” column).\ More details about the prediction accuracy abs_R has been provided in *[Prediction_accuracy_R_data.csv](https://www.dropbox.com/scl/fi/fh0adwhqc7higr26jaan1/Prediction_accuracy_R_data.csv?rlkey=krsxawc05jups57n40ebpmw2y&st=u2ukx6bu&dl=1)*.
