####
args = commandArgs(trailingOnly=TRUE)
k <- args[1]
k=as.numeric(k)

library(data.table)
res=data.frame(PF=NA,eQTL=NA,eGene=NA,eQTL_raw=NA,eGene_raw=NA)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

#freq=c()

#for(i in 1:22){

#  tmp=fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/reverse_vQTL/OneK1K_SNPs_freq_chr",i,".txt"),header=T)
#  tmp=as.data.frame(tmp)

#  freq=rbind(freq,tmp)

#  print(i)
#}

#low=freq[freq$freq<0.05,]

# Pi for each gene
# gene=read.table(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/mean_variance_relationship/linear_reg/",new_names[k],"_cells_mean_var_cor_and_e_normal_test.txt"),header=T)

# pi=gene[gene$Pi_wi_var>0.9,]

for(i in 0:50){

a=fread(paste0("./",new_names[k],"/",i,"/",new_names[k],"_eQTL_all.txt"),header=T)
a=as.data.frame(a)

res[i+1,4]=nrow(a)
res[i+1,5]=length(unique(a$gene))

# a=a[a$LFDR<0.05,]
# a2=a
# a2=a[!a$snps %in% low$snps,]
# a3=a2[!a2$gene %in% pi$gene,]

res[i+1,1]=i
res[i+1,2]=nrow(a)
res[i+1,3]=length(unique(a$gene))
print(i)
}

write.table(res,paste0(new_names[k],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_all_genes_PF1_50.txt"),row.names=F,col.names=T,quote=F,sep="\t")


####
