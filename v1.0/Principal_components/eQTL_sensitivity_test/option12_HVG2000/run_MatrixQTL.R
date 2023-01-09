##############################################################################
# Script information
# Title: Sensitivity test for eQTL analysis in OneK1K
# Author: Angli Xue
# Date: 2022-09-24
# Description: This R script was written to run the MatrixeQTL analysis with 0-50 PEER factors for 14 cell types in OneK1K
##############################################################################

print(Sys.time())
start=Sys.time()

args = commandArgs(trailingOnly=TRUE)

chrNumber <- args[1]
PF_num <- args[2]
ct_name <- args[3]

ct_name=as.character(ct_name)
PF_num=as.numeric(PF_num)


system(paste0("mkdir -p /share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/PF_test/PCA/HVG2000/",ct_name,"/"))
system(paste0("mkdir -p /share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/PF_test/PCA/HVG2000/",ct_name,"/",PF_num,"/"))
setwd(paste0("/share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/PF_test/PCA/HVG2000/",ct_name,"/",PF_num))

library("MatrixEQTL")
library(qvalue)

# Cell types
cell <- c("TCL1A+FCER2+Bcell", "TCL1A-FCER2-Bcell", "CD4+KLRB1-Tcell", "CD4+KLRB1+Tcell",
    "CD4+SOX4+Tcell", "CD8+LTB+Tcell", "CD8+GNLY+NKG7+Tcell", "CD8+S100B+Tcell",
    "Dendriticcell", "MonocyteCD14+", "MonocyteFCGR3A+", "XCL1+NK", "XCL1-NK", "IgJ+Bcell",
    "Erythrocytes","Platelets")
new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
    "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
    "Erythrocytes","Platelets")

new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff",
               "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma",
               "Erythrocytes","Platelets")

base.dir = paste0("/share/ScratchGeneral/seyyaz/onek1k/cell_specific_eQTL_analysis_October19/",cell[which(new_names==ct_name)],"/")

base.dir2 = paste0("/directflow/SCCGGroupShare/projects/SeyhanYazar/onek1k/science_revision_Dec20/peer_factors/",new_names2[ct_name],"/")

# Then we set the parameters such as selected linear model and names of genotype and expression data files.
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
# Genotype file name
SNP_file_name = paste0(base.dir,'step1/SNPs_chr',chrNumber,'.txt', collapse='');
snps_location_file_name = paste0(base.dir,'step1/snpsloc_chr',chrNumber,'.txt', collapse='');

# Gene expression variance file name
# system(paste0("mkdir -p /share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/input_files/",ct_name,"/"))
exp.dir=paste0("/share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/mean_mx_input_files/",ct_name,"/")
expression_file_name = paste(exp.dir, "expression_chr",chrNumber,".txt", sep="");
gene_location_file_name = paste0(base.dir, 'step1/geneloc_chr',chrNumber,'.txt', collapse='');

covariates_file_name = paste0("/share/ScratchGeneral/angxue/proj/vQTL/MatrixQTL/round2_more_covar/covariates/no_outliers_PF50_HVG2000/PCA/",new_names[which(new_names==ct_name)],"_PC",PF_num,".txt", collapse='');

# Output file name
output_file_name_cis = paste0('eQTL_chr',chrNumber,'.txt');
#output_file_name_tra = tempfile();


# threshold
# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
#pvOutputThreshold_tra = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric()

# Distance for local gene-SNP pairs
cisDist = 1e6;

# Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Match column ID between snps and gene
if(length(gene$columnNames) < length(snps$columnNames)){
index=match(gene$columnNames,snps$columnNames)
snps$ColumnSubsample(index)
}

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

#
## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

# main function
me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
#output_file_name     = output_file_name_tra,
#pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
unlink(output_file_name_cis);


## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');

print(head(me$cis$eqtls))
print(paste0(nrow(me$cis$eqtls)," SNP-gene pairs tested."))

#cat('Detected distant eQTLs:', '\n');
print("Calculating the local FDR...")
qobj = qvalue(p = me$cis$eqtls$pvalue)
summary(qobj)

me$cis$eqtls$qvalue = qobj$qvalues
me$cis$eqtls$LFDR = qobj$lfdr

## Plot the Q-Q plot of local and distant p-values
sig=me$cis$eqtls
sig=sig[sig$LFDR<0.05,]
sig=sig[order(sig$snps),]

print(paste0("Identified number of eGenes: ",length(unique(sig$gene))))
#plot(me)
write.table(sig,file=paste0(ct_name,'_eQTL_chr',chrNumber,'_FDR_005.txt'),quote=F,row.names=F,sep="\t")

# significant_cis = me$cis$eqtls
# write.table(significant_cis,file=paste0(ct_name,'_vQTL_chr',chrNumber,'_FDR_0.txt'),quote=F,row.names=F,sep="\t")

print(paste0("chr",chrNumber," for ",ct_name," cells finished!"))

print(Sys.time())
end = Sys.time()
z = as.numeric(end - start, units = "mins")
print(paste0("Running time for whole script is ", z, " mins"))


#### END
