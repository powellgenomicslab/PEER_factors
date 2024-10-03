####
#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(paste0("Read in argument, ",args))

library(data.table)

cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell",
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell",
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")
new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

ct_name = as.character(new_names[as.numeric(args[1])])

system(paste0("mkdir -p /share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/mean_mx_input_files/",ct_name,"/"))
setwd(paste0("/share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/mean_mx_input_files/",ct_name,"/"))

mean_mx=fread(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/",ct_name,"_cells_mean_mx.txt"),header=T)
mean_mx=as.data.frame(mean_mx)

# Exclude genes with pi0>0.9
pi0=rowSums(mean_mx[,-1]==0,na.rm=T)/(ncol(mean_mx)-length(which(is.na(colSums(mean_mx)))))

mean_mx=mean_mx[which(pi0<=0.9),]

# log+1 transform
tmp=colSums(mean_mx)
excol=which(is.na(tmp))
mean_mx[,-excol]=log1p(mean_mx[,-excol])

# Z-score standardization
mean_mx=t(scale(t(mean_mx)))

gene=read.table(paste0("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/",ct_name,"_cells_gene_list.txt"),header=T)

mean_mx=cbind(gene[which(pi0<=0.9),],mean_mx)
# mean_mx=cbind(gene,mean_mx)
colnames(mean_mx)[1]="id"
mean_mx=as.data.frame(mean_mx)

## read Seyhan's format
base_dir = paste0("/share/ScratchGeneral/seyyaz/onek1k/cell_specific_eQTL_analysis_October19/",cell[which(new_names==ct_name)],"/step1/")

for (i in 1:22){
exp=fread(paste0(base_dir,"expression_chr",i,".txt"),header=T)
exp=as.data.frame(exp)

# mean_mx=mean_mx[,colnames(exp)]
all(colnames(mean_mx) %in% colnames(exp))

# align the mean_mx and exp
extra=colnames(exp)[!colnames(exp) %in% colnames(mean_mx)]
v2=matrix(NA,nrow=nrow(mean_mx),ncol=length(extra))
v2=as.data.frame(v2)
colnames(v2)=extra
mean_mx=cbind(mean_mx,v2)
mean_mx=mean_mx[,match(colnames(exp),colnames(mean_mx))]
all(colnames(exp)==colnames(mean_mx))

g=Reduce(intersect,list(exp$id,mean_mx$id))

# Subset gene list
mean_mx_new=mean_mx[match(g,mean_mx$id),]

print(paste0(nrow(mean_mx_new)," genes matched."))
write.table(mean_mx_new,paste0("expression_chr",i,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
print(paste0(ct_name," cells, chr ",i," format conversion finished!"))
}

#### END
