##############################################################################
# Script information
# Title: Downsampling of iPSCs data in Neavin et al.2021.Genome Biology.
# Author: Angli Xue
# Date: 2022-01-25
# Description: This R script was written to run a job for downsampling analysis
##############################################################################

args = commandArgs(trailingOnly=TRUE)

rep <- args[1]
rep=as.numeric(rep)

# Import libraries
library(peer)
library(data.table)
library(RNOmni)

set.seed(2022*rep)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff",
               "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma",
               "Erythrocytes","Platelets")


# Get expression file --------------------------------
expr=read.table("Expression4PEER_cluster1.tsv",header=T)
dim(expr)
rnames = as.character(rownames(expr))

# Change N to 40 and 50 for different downsampling size
N=31

id=sample(1:nrow(expr),N)
expr2=expr[id,]

#Set PEER paramaters based on the instructions from PEER package website

model = PEER()

PEER_setPhenoMean(model, as.matrix(expr2))

# PEER_setCovariates(model, as.matrix(covs2))

PEER_setNmax_iterations(model,1000)

dim(PEER_getPhenoMean(model))

# PEER_setAdd_mean(model, TRUE)

# PEER_setNmax_iterations(model, 100)

pf_num=10 # Number of PEER factors

PEER_setNk(model,pf_num) # Set to generate 10 PEER factors

PEER_getNk(model)

PEER_update(model)

# saveRDS(model,"pf50.rds")

# Visualiza the variance per PEER factor
alpha = PEER_getAlpha(model)
# alpha=alpha[1:length(alpha),]
write.table(1/alpha,paste0("cluster1_weights_PF",pf_num,"_rep",rep,".txt"),row.names=F,col.names="Variance",quote=F)

# pdf(file=paste0(ct_name,"_plot_precision_HVG2000_PF",pf_num,"_RINT_scaled.pdf"),width=8,height=7)
# plot(1/alpha,xlab="Factors",ylab="Variance of factor weights", type="b", col="blue", lwd=4, xaxp=c(1,length(alpha), length(alpha)-1),main=ct_name)
# dev.off()

# Calculate and save the PEER factors
factors = PEER_getX(model)
dim(factors)
factors_df = data.frame(factors)
factors_df$sampleid = rnames[id]
# write.table(factors_df,"test_addMean.txt",col.names=T,row.names=F,quote=F)
factors_df = factors_df[c((pf_num+1),1:(pf_num))]
colnames(factors_df) = c(c("sampleid"),paste0("PF",1:pf_num))

write.table(factors_df,paste0("Expression4PEER_cluster1_PF10_rep",rep,".txt"),col.names=T,row.names=F,quote=F)

corr=cor(factors_df[,2:11])
message(formatC(mean(abs(corr[lower.tri(corr)])),digits=4))




