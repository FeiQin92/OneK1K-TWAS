Files required for runnning FUSION-TWAS

1. Weights matrix (including 14 cell types in our data)
 
----Cell type: "CD4_ET", "NK", "CD4_NC", "CD8_S100B", "CD8_ET", "B_IN", "CD8_NC", "B_Mem", "NK_R", "Mono_NC", "Mono_C", "DC", "Plasma", "CD4_SOX4".
----To identify more susceptibility genes through TWAS, we implemented three distinct elastic net-based approaches to build gene expression prediction models for each cell type: 1) the traditional 'targetC'ù model, which predicts expression using data only from the targeted cell type; 2) the 'S+targetC'ù model, following Thompson et al., which incorporates both shared and cell type-specific expression components for the targeted cell type; and 3) a novel 'S+allC'ù model proposed by us, which integrates shared and cell type-specific expression components across all cell types to further improve prediction performance in the targeted cell type.
----In the weights files, which have been mangaed according to the format requirements of FUSION, weights for each of three models (targetC, S+targetC, S+allC) above have been provided. The SNP id is express as "chr_pos", i.e., 1_918384.



2. GWAS summmary statistics, organized in a file including columns of
#SNP: SNP id with format of "chr_pos". We utlized postion from hg19 reference.
#A1: Alternative allele;
#A2: Reference alllel; 
#Z: GWAS Z statistics  

Example:
#SNP     A1      A2      Z
#1_918384        T       G       0.307692307692308
#1_918573        G       A       0.0879120879120879
#1_962606        A       G       1.96666666666667 

#Below are my example R script to generate this GWAS summary file
#BC_GWAS <- data.frame(SNP=paste0(GWAS2$Chr, "_", GWAS2$Pos),
#                        A1=GWAS2$ALT, 
#                        A2=GWAS2$Ref, 
#                        Z=GWAS2$beta/GWAS2$se)
#write.table(BC_GWAS, file="..../BC_fusion.txt", quote=F, row.names=F, col.names=T, sep="\t")



3. refernece genome
We utlized hg19 EUR from 1000 Genome. You can use the one you want.



4. Analysis code for FUSION_TWAS for each of three prediction models were also provided.
----FUSION_targetC.R
----FUSION_targetC.swarm
----FUSION_S_targetC.R
----FUSION_S_targetC.swarm
----FUSION_S_allC.R
----FUSION_S_alltargetC.swarm
The only difference for three R files is the "force_model" in the function.
In each .swarm file, the 2th column is the index for chromosome, the 3th column is the index for cell type.  