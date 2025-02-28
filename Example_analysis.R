library(bigstatsr)
library(data.table)
library(caret)
library(glmnet)

######################## Below are scripts to generating weights matrix required for FUSION TWAS #################
##################################################################################################################


setwd(paste0("C:/Users/qinf2/OneDrive - National Institutes of Health/Kai Yu/eQTL/TWAS/Code/Github OneK1K_TWAS/code"))

source("1_hom_het_components_fit.R")
source("2_Model_fit_using_hom_het_components.R")
source("3_Weights_Pre.R")
source("4_Gen_weights.R")

setwd(paste0("C:/Users/qinf2/OneDrive - National Institutes of Health/Kai Yu/eQTL/TWAS/Code/Github OneK1K_TWAS/Data"))


# type="logCPM"
# ID_for_swarm=args[1] #Gene list file 000:498

############################ Fiting hom and het components for each cell type #############################
###########################################################################################################
###########################################################################################################
# genelist <- read.table(paste0("/gpfs/gsfs12/users/qinf2/OneK1K/TWAS/Gelist_50_all", ID_for_swarm))
# for (GID in 1:nrow(genelist)){
#   genename <- genelist$V1[GID]
  
  genename <- "ENSG00000000457"
  X1 <- paste0("Exp/", genename, "/snpexp")
  Y1 <- paste0("Exp/", genename, "/y/")
  cov1 <- "Cov/"
  Out1 <- "results/"
  SNP1 <- paste0("Exp/", genename, "/snppos.txt")
  Pre1 <- genename
  
  # print(paste0(GID, "th Gene: ", genename))
  
  # 
  # X=NULL; 
  # Ys=NULL; 
  # covs=NULL; 
  # X_file=X1;  # Genotype matrix, individuals are row names, no column names.
  # Y_file_dir=Y1; # Direction containing only the expression files, individuals are row names, no column names.
  # cov_file_dir=cov1; # Cov file direction, should share the same prefix as the expression files
  # verbose=T;
  # out_dir=Out1;  # output CV predictors and other statistics
  # snps=SNP1;  # A tab-delimited txt file containing snp information. six columns.
  # signif=0.1;
  # gene_name=Pre1; # Add this prefix to all the saved results files
  # content_alpha=0.5;  # Regularization constant, 1 means LASSO regression.
  # tissue_alpha=0.5;  # Regularization constant for the tissue by tissue approach.
  # num_folds=10;  # Number of folds for cross-validation
  # scale_exp=F;  # Center and scale the gene expression
  # seed=9000
  
  hom_het_model(X_file=X1,  # Genotype matrix, individuals are row names, no column names.
                Y_file_dir=Y1, # Direction containing only the expression files, individuals are row names, no column names.
                cov_file_dir=cov1, # Cov file direction, should share the same prefix as the expression files
                out_dir=Out1,  # output CV predictors and other statistics
                snps=SNP1,  # A tab-delimited txt file containing snp information. six columns.
                gene_name=Pre1) # Add this prefix to all the saved results files
# }	



# SNP1 <- read.table(paste0("Exp/", genename, "/snppos.txt"))
# SNP1$V2 <- paste0(SNP1$V1, "_" ,SNP1$V4)
# a <- table( SNP1$V2)
# dup_names <- names(a)[a>1]
# SNP2  <- SNP1[!SNP1$V2 %in% dup_names, ]
# write.table(SNP2, file=paste0("Exp/", genename, "/snppos.txt"), row.names = F, col.names = F, sep = "\t", quote=F)
# 
# X1 <- read.table(paste0("Exp/", genename, "/snpexp"))
# X2 <- X1[, c(TRUE, !SNP1$V2 %in% dup_names)]
# write.table(X2, file=paste0("Exp/", genename, "/snpexp"), row.names = F, col.names = F, sep = "\t", quote=F)

################################### Second step ######################################
######  Fitting models using hom and het components generated in the first step ######
######################################################################################  
Model_fit(Gname=genename)

  
######################################################################################
###### Combing coefficients from above two steps and prepare for weights file   ######
######################################################################################    
Weights_pre(G_name=genename, in_dir="results/")






######################################################################################
############# Generating weights files required for FUSION TWAS ######################
######################################################################################   
# OneK1K_ldprunekeep <- data.table::fread("/gpfs/gsfs12/users/qinf2/OneK1K/TWAS/logCPM/model/GWAS/OneK1K_ldprunekeep.txt", header=F)
# snp_pos = data.table::fread('/gpfs/gsfs12/users/qinf2/OneK1K/TWAS/snp_pos.csv', header = TRUE, sep = ",")  # Replace with your SNP data file
# snp_pos <- as.data.frame(snp_pos)
# snp_pos$SNP <- paste0(snp_pos$chromosome, "_", snp_pos$Position)
# rownames(snp_pos) <- paste0(snp_pos$chromosome, "_", snp_pos$Position, "_", snp_pos$A1, "_", snp_pos$A2)
# snp_pos$Chr <- snp_pos$chromosome
# snp_pos$Pos <- snp_pos$Position
# snp_pos_sub <- snp_pos[OneK1K_ldprunekeep$V1,]
# rownames(snp_pos_sub) <- snp_pos_sub$SNP
Weights_gen(G_name=genename, in_dir="results/", out_dir="weights/")
    