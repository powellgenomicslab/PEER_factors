####
library(ggplot2)
library(tidyr)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
               "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
               "Erythrocytes","Platelets")

# new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R","Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col=c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")


a=c()

for(i in 1:14){
  #b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_all_genes_PF1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  b$ct_col=ct_col[i]
  a=rbind(a,b)
}

old=a

a=c()

for(i in 1:14){
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_HVG2000_PF1_50.txt"),header=T,sep="\t")
  # b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/new/17APR2022/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  b$ct_col=ct_col[i]
  a=rbind(a,b)
}

new=a

old$Group="All_genes"
new$Group="HVG2000"

old=rbind(old,new)

old$ct = factor(old$ct, levels = new_names[1:14])

##
## Calculate how much signals can fitting PFs gain
st1=data.frame(ct=new_names[1:14],N=NA,All_genes=NA,HVG2000=NA,nr_max1=NA,nr_max2=NA)

for(i in 1:14){
  # All genes
  tmp = old[old$ct==new_names[i] & old$Group=="All_genes",]
  base = tmp$eGene[1]
  p.incr = max(tmp$eGene[-1])
  per = (p.incr - base) / base
  st1$All_genes[i] = formatC(per,format="f",digits=4)
  st1$nr_max1[i]=which.max(tmp$eGene[-1])
  
  # HVG2000
  tmp = old[old$ct==new_names[i] & old$Group=="HVG2000",]
  base = tmp$eGene[1]
  p.incr = max(tmp$eGene[-1])
  per = (p.incr - base) / base
  st1$HVG2000[i] = formatC(per,format="f",digits=4)
  st1$nr_max2[i]=which.max(tmp$eGene[-1])
}


ct_sam=read.table("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Table/nr_ind_by_QC_min_cells.txt",header=T)
ct_1=read.table("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Table/nr_ind_by_QC_min_cells_v2.txt",header=F,sep=" ")
ct_1=ct_1[,-1,drop=F]

ct_1=separate(ct_1,V2,into=paste0("V",1:10), sep=" ")
ct_1=ct_1[,c(6,9)]
colnames(ct_1)=c("cell_type","ct1")
ct_1=ct_1[match(ct_sam$cell_type,ct_1$cell_type),]
ct_sam$ct1=ct_1$ct1
ct_sam=ct_sam[,c(1,5,2:4)]

all(st1$ct==ct_sam$cell_type)
st1$N=ct_sam$ct1

write.table(st1,"~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Table/Table_S1_eGene_power_gain.txt",row.names = F,col.names = T,quote = F, sep="\t")


