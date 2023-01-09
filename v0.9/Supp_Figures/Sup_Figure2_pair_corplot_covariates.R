#### 
library(psych)


new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
               "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
               "Erythrocytes","Platelets")
new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff",
                "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma",
                "Erythrocytes","Platelets")

new_names3 <- c(expression(B[IN]), expression(B[Mem]), expression(CD4[NC]), expression(CD4[ET]), expression(CD4[SOX4]), 
                expression(CD8[NC]), expression(CD8[ET]), expression(CD8[S100B]), "DC", 
                expression(Mono[C]), expression(Mono[NC]), expression(NK[R]), "NK", "Plasma",
                "Erythrocytes","Platelets")
k=3

## Figure 1A 
old=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_1/",new_names[k],"_mx_pf50.txt"),header=T,check.names = F)

old[,-1] = as.data.frame(sapply(old[,-1], as.numeric))
colnames(old)=c("ID","Sex",paste0("PC",1:6),"Age",paste0("PF",1:50))

pdf(file="~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure2_pair_corrplot_covariates_v2.pdf",width = 9,height = 7)

pairs.panels(old[,2:9], 
             hist.col = viridis::viridis(10)[4],
             smooth = F,
             stars=T,
             ellipses = F, # show correlation ellipses
             cex = 0.1,
             cex.labels = 0.9,
             cex.cor = 1, xaxt = 'n', yaxt = 'n', alpha=0.02)

dev.off()







####