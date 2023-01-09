####
library(data.table)
cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell",
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell",
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")
#

new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff",
               "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma",
               "Erythrocytes","Platelets")

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")


base.dir = paste0("/share/ScratchGeneral/angxue/proj/vQTL/PEER_factors/scale_center_test/no_outliers/option_1/")

for(ct_name in 1:14){
a=read.table(paste0(base.dir,new_names[ct_name],"_mx_pf50.txt"),header=T)

a2=t(a)
a2=as.data.frame(a2)

colnames(a2)=a2[1,]
a2=a2[-1,]

a2=cbind(rownames(a2),a2)
colnames(a2)[1]="id"

base.dir2 = paste0("/share/ScratchGeneral/seyyaz/onek1k/cell_specific_eQTL_analysis_October19/",cell[ct_name],"/")
snp=fread(paste0(base.dir2,'step1/SNPs_chr22.txt', collapse=''))
snp=as.data.frame(snp)
a2[colnames(snp)[which(!colnames(snp) %in% colnames(a2))]] <- NA
a2=a2[,match(colnames(snp),colnames(a2))]

for(i in 0:50){

tmp=a2[1:(8+i),]

write.table(tmp,paste0(new_names[ct_name],"_covar_peer_factors_PF",i,".txt"),row.names=F,col.names=T,quote=F,sep="\t")

}

print(new_names[ct_name])

}

####
