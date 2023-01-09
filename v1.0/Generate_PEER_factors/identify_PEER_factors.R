##############################################################################
# Script information
# Title: Identify PEER factors in CD4 NC cells
# Author: Angli Xue
# Date: 2022-01-25
# Description: This R script was written to run a job for PEER factor identification
# This script was modified based on the script from Seyhan Yazar
##############################################################################

print(Sys.time())
start=Sys.time()

args = commandArgs(trailingOnly=TRUE)

ct_name <- args[1]
ct_name=as.numeric(ct_name)

opt <- args[2]
opt = as.numeric(opt)

system(paste0("mkdir -p /share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/option_",opt))
setwd(paste0("/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/option_",opt))

# Import libraries
library(peer)
library(data.table)
library(RNOmni)
library(moments)
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
param=readRDS("/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/param_options.RDS")
# expr <- pseudobulk_scaling(expr=expr,pi0=0.9,log1p=T,scale=T,RINT=F,HVG=2000)
Sys.sleep(runif(1))
expr2 <- pseudobulk_scaling(expr=expr,
	pi0=param[[opt]]$pi0,
	log1p=param[[opt]]$log1p,
	scale=param[[opt]]$scale,
	RINT=param[[opt]]$RINT,
	HVG=param[[opt]]$HVG)

sk=apply(expr2,2,function(x){skewness(x)})
sk=data.frame(gene=colnames(expr2),skewness=sk)
write.table(sk,paste0(ct_name,"_mx_skewness.txt"),col.names=T,row.names=F,quote=F)


#Set PEER paramaters based on the instructions from PEER package website

model = PEER()

PEER_setPhenoMean(model, as.matrix(expr2))

PEER_setCovariates(model, as.matrix(covs2))

PEER_setNmax_iterations(model,2000)

dim(PEER_getPhenoMean(model))

# PEER_setAdd_mean(model, TRUE)

# PEER_setNmax_iterations(model, 100)

pf_num=50 # Number of PEER factors

PEER_setNk(model,pf_num) # Set to generate 10 PEER factors

PEER_getNk(model)

print("Start the PEER estimation...")
print(Sys.time())
tick1=Sys.time()

PEER_update(model)

print(Sys.time())
tick2=Sys.time()
z=as.numeric(tick2-tick1,units="mins")
print(paste0("Running time for PEER estimation is ",z," mins"))

# saveRDS(model,"pf50.rds")

# Visualiza the variance per PEER factor
alpha = PEER_getAlpha(model)
alpha=alpha[9:length(alpha),]
write.table(1/alpha,paste0(ct_name,"_variance_pf",pf_num,".txt"),row.names=F,col.names="Variance",quote=F)

pdf(file=paste0(ct_name,"_plot_precision_pf",pf_num,".pdf"),width=8,height=7)
plot(1/alpha,xlab="Factors",ylab="Weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1),main=ct_name)
dev.off()

# Calculate and save the PEER factors
factors = PEER_getX(model)
dim(factors)
factors_df = data.frame(factors)
factors_df$sampleid = rnames
# write.table(factors_df,"test_addMean.txt",col.names=T,row.names=F,quote=F)
factors_df = factors_df[c((pf_num+9),1:(pf_num+8))]
colnames(factors_df) = c(c("sampleid", "sex", "pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "age"),paste0("pf",1:pf_num))

# corr=cor(factors_df[,10:ncol(factors_df)])

write.table(factors_df,paste0(ct_name,"_mx_pf",pf_num,".txt"),col.names=T,row.names=F,quote=F)

# residuals = PEER_getResiduals(model)
# dim(residuals)
# residuals_df = data.frame(residuals)
# colnames(residuals_df) = colnames(expr2)
# residuals_df = cbind(sampleid=rnames, residuals_df)
# write.table(residuals_df,paste0(ct_name,"_residuals_pf",pf_num,".txt"),col.names=T,row.names=F,quote=F)

print(Sys.time())
end=Sys.time()
z=as.numeric(end-start,units="mins")
print(paste0("Running time for whole script is ",z," mins"))



