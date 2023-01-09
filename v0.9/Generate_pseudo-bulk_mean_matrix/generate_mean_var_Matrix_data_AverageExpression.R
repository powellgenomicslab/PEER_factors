#### OneK1K data set ####
# CD4_NC_raw_counts.RDS  ---- original counts
# CD4_NC_rna_counts.RDS  ---- log normalised counts
# CD4_NC_sct_counts.RDS  ----  corrected counts for sequencing depth
# meta_data.RDS ---- this file has the individual information. Ie which barcode belongs to which individual.

print(Sys.time())
start=Sys.time()

#! /usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(paste0("Read in argument, ",args))
library(MASS)
setwd("/share/ScratchGeneral/angxue/proj/vQTL/mean_var_matrix/")

# Cell type names
new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell",
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell",
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")

ct_name = as.character(new_names[as.numeric(args[1])])

print(ct_name)

## Read the data ##
info=readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/",ct_name,"_meta_data.RDS"))
dim(info)
## Number of cells and individuals
print(paste("This data set includes",dim(info)[1],ct_name,"cells from",nrow(table(info$individual)),"individuals"))

print("Number of cells per individual: ")
print(summary(as.numeric(table(as.character(info$individual)))))

## Read sctransformed counts 
print("Start to read in sct counts...")
sct=readRDS(paste0("/directflow/SCCGGroupShare/projects/angxue/data/onek1k/",ct_name,"_sct_counts.RDS"))
# Check data structure
str(sct)

## Load package
library(dplyr)
library(Seurat)

# Convert sct counts to Seurat object
ok1k <- CreateSeuratObject(counts = sct, assay = "SCT", project = "onek1k", min.cells = 0, min.features = 0)
ok1k

# Match Seurat object with previously saved meta.data
# Number might be inconsistent due to filters (e.g., min.cells, min.features)
info=info[info$UMI %in% rownames(ok1k@meta.data),]

# Add back individual ID from meta.data 
ok1k@meta.data$individual=info$individual

# Extract a list of individual ID and their number of cells
a=unlist(table(info$individual))

## Check mean-variance relationship within an individual
# To ensure the estimated within-ind variance is accurate, we excluded those individuals with < 3 cells (this threshold needs better check)
# keep_list=names(which(a>=3))
# ok1k=subset(x = ok1k, subset = individual %in% keep_list)
# print(paste0(length(keep_list)," individuals remain after excluding those with < 3 cells."))
# Update the meta.data (info)
# info=info[info$UMI %in% rownames(ok1k@meta.data),]
print(summary(as.numeric(table(as.character(info$individual)))))

## Main function to calculate the pseudobulk mean expression
Idents(object=ok1k ) <- 'individual'
length(levels(x=ok1k))

person.averages <- AverageExpression(object= ok1k, assays = "SCT", features = NULL,
  return.seurat = FALSE, group.by = "individual", slot = "counts",
  use.scale = FALSE, use.counts = FALSE, verbose = TRUE)

person.averages=as.data.frame(person.averages[[1]])

gene=rownames(person.averages)

# print(paste0(length(gene)," genes included in the main analysis"))
write.table(as.data.frame(gene),paste0(ct_name,"_cells_gene_list.txt"),row.names=F,col.names=T,quote=F)

# Aggregate Pseudobulk
person.sums <- AggregateExpression(object= ok1k, assays = "SCT", features = NULL,
  return.seurat = FALSE, group.by = "individual", slot = "counts",
  use.scale = FALSE, use.counts = FALSE, verbose = TRUE)

person.sums=as.data.frame(person.sums[[1]])

gene=rownames(person.sums)


# Match the column with covarites file
co=read.table(paste0("/share/ScratchGeneral/seyyaz/onek1k/cell_specific_eQTL_analysis_October19/",cell[which(new_names==ct_name)],"/step1/covariates_chr1.txt"),header=T,check.names=F)
all(colnames(person.averages) %in% colnames(co)[-1])

extra=colnames(co)[!colnames(co) %in% colnames(person.averages)]
if(length(extra)>1){

extra=extra[-1]
tmp=matrix(NA,nrow=nrow(person.averages),ncol=length(extra))
tmp=as.data.frame(tmp)
colnames(tmp)=extra

person.averages=cbind(person.averages,tmp)
person.sums=cbind(person.sums,tmp)

print(paste0(ncol(co)-1," individuals matched with covariate files."))
all(colnames(person.averages) %in% colnames(co)[-1])
all(colnames(co)[-1] %in% colnames(person.averages))
person.averages=person.averages[,colnames(co)[-1]]
} else{
person.averages=person.averages[,colnames(co)[-1]]
person.sums=person.sums[,colnames(co)[-1]]
}

if(all(colnames(person.averages)==colnames(co)[-1])){print("Column match checked!")}

## Save results
write.table(person.averages,paste0(ct_name,"_cells_mean_mx.txt"),row.names=F,col.names=T,quote=F)

write.table(person.sums,paste0(ct_name,"_cells_sum_mx.txt"),row.names=F,col.names=T,quote=F)

print("Script ends.")
print(Sys.time())
end=Sys.time()
print(difftime(end, start, units = "mins")  )

#### END ####
