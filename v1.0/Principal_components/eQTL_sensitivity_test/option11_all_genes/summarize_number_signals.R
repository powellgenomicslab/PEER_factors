####
args = commandArgs(trailingOnly=TRUE)
k <- args[1]
k = as.numeric(k)

library(data.table)
res=data.frame(PF = NA, eQTL = NA, eGene = NA, eGene_tested = NA)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

for(i in 0:50){

a=fread(paste0("./",new_names[k],"/",i,"/",new_names[k],"_eQTL_all.txt"),header=T)
a=as.data.frame(a)

res[i+1,1]=i
res[i+1,2]=nrow(a)
res[i+1,3]=length(unique(a$gene))
# print(i)
}

# Calculated the total number of genes tested 
nr_gene = system(paste0("wc -l /share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/mean_mx_input_files/",new_names[k],"/*.txt | awk '/ total/{ print $1 }'"), intern = T)
nr_gene = as.numeric(nr_gene) - 22

res$eGene_tested = nr_gene

write.table(res, paste0(new_names[k],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_all_genes_PC1_50.txt"), row.names = F, col.names = T, quote = F, sep = "\t")

print(new_names[k])

####
