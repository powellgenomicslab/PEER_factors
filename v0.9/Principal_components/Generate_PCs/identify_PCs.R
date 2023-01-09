##############################################################################
# Script information
# Title: Identify PCs in CD4 NC cells
# Author: Angli Xue
# Date: 2022-09-25
# Description: This R script was written to run a job for PCA anlaysis
##############################################################################

print(Sys.time())
start=Sys.time()

args = commandArgs(trailingOnly=TRUE)

ct_name <- args[1]
ct_name=as.numeric(ct_name)

opt <- args[2]
opt = as.numeric(opt)

system(paste0("mkdir -p /share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/PCA/option_",opt))
setwd(paste0("/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/PCA/option_",opt))

# Import libraries
library(PCAForQTL)
library(data.table)
library(RNOmni)
# library(moments)
source("/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/pseudobulk_scaling.R")

set.seed(2022)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

# new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff",
#                "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma",
#                "Erythrocytes","Platelets")
ct_name=new_names[ct_name]
print(paste0("Cell type: ",ct_name))

# Get expression file --------------------------------
expr=fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/",ct_name,"_cells_mean_mx.txt"),header=T)
expr=as.data.frame(expr)

dim(expr)
gene=read.table(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/",ct_name,"_cells_gene_list.txt"),header=T)
expr=as.data.frame(t(expr))
colnames(expr)=gene$gene

# Remove individuals with NA
na_idx=which(is.na(rowMeans(expr)))
if(length(na_idx) > 0){
expr = expr[-na_idx,]
}

dim(expr)

expr[1:5,1:5]
rnames=rownames(expr)
cnames = colnames(expr)

# Get covariate files ---------------------------------
covs = read.csv('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/covariates.tsv', sep='\t')
age = read.csv('/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/cell_specific_eQTL_analysis_April20/age_covariate.tsv', sep='\t')
colnames(age)[1] <- "sampleid"

# Merge covarite and age file
covs = merge(covs, age, by="sampleid")
dim(covs)

# Check all the individuals exist in both files
table(rnames %in% covs$sampleid)

# ***IMPORTANT*** Match the order of individuals
covs = covs[match(rnames, covs$sampleid),]
covs2 = covs[,-1]

# Transform and scale the matrix
param = readRDS("/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/param_options.RDS")
# expr <- pseudobulk_scaling(expr=expr,pi0=0.9,log1p=T,scale=T,RINT=F,HVG=2000)
Sys.sleep(runif(1))
expr2 <- pseudobulk_scaling(expr=expr,
	pi0=param[[opt]]$pi0,
	log1p=param[[opt]]$log1p,
	scale=param[[opt]]$scale,
	RINT=param[[opt]]$RINT,
	HVG=param[[opt]]$HVG)


print("Start the PCA estimation...")
print(Sys.time())
tick1 = Sys.time()

prcompResult <- prcomp(expr2, center = F, scale = F) #This should take less than a minute.
PCs <- prcompResult$x

print(Sys.time())
tick2 = Sys.time()
z = as.numeric(tick2-tick1,units="mins")
print(paste0("Running time for PCA estimation is ",z," mins"))

# Add known covariates back
PC_num = 50
factors_df = cbind(sampleid = rownames(PCs), covs2, PCs[,1:PC_num])

colnames(factors_df) = c(c("sampleid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age"), paste0("PCA", 1:PC_num))

write.table(factors_df, paste0(ct_name, "_mx_pca", PC_num, ".txt"), col.names = T, row.names = F, quote = F)

print(Sys.time())
end = Sys.time()
z = as.numeric(end - start, units = "mins")
print(paste0("Running time for whole script is ", z, " mins"))



