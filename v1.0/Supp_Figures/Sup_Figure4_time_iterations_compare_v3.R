#### Sup Table or Figure for time and iterations
library(ggplot2)
library(tidyr)
library(readr)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
               "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
               "Erythrocytes","Platelets")

new_names2 <- c("BimmNaive", "Bmem", "CD4all", "CD4effCM", "CD4TGFbStim", "CD8all", "CD8eff",
                "CD8unknown", "DC", "MonoC", "MonoNC", "NKact", "NKmat", "Plasma",
                "Erythrocytes","Platelets")

new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R","Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col=c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")


setwd("/Users/uqaxue/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/")

# case=c("pf10","pf10_pi1_excluded","pf10_pi09_excluded","pf10_log1p_scaled_pi1_excluded","pf10_log1p_scaled_pi09_excluded","pf10_scaled_pi1_excluded","pf10_scaled_pi09_excluded")

res=c()

for(i in 6:13){
  
  t=read.table(paste0("./option_",i,"/time_option_",i,".txt"),header=F)
  t=t[grep("PEER",t$V2),]
  t$V1=gsub(".*[.]([^.]+)[:].*", "\\1", t$V1)
  t$V1=as.numeric(t$V1)
  t=t[order(t$V1,decreasing = F),]
  t$V2=parse_number(t$V2)
  t$V1=new_names[1:14]
  t$option=i
  
  res=rbind(res,t)
  
}
colnames(res)=c("cell_type","time","option")

res$ct_col = ct_col
res$ct_col = factor(res$ct_col, levels = ct_col)


#### Iterations
res2=c()

for(i in 6:13){
  t=read.table(paste0("./option_",i,"/iterations_",i,".txt"),header=F)
  t$V1=gsub(".*[.]([^.]+)[:].*", "\\1", t$V1)
  t$V1=as.numeric(t$V1)
  t=t[,c(1,4)]
  if(i==1){t=rbind(t,c(9,1000))}
  t=t[order(t$V1,decreasing = F),]
  t$V1=new_names[1:14]
  t$option=i
  
  res2=rbind(res2,t)
}

colnames(res2)=c("cell_type","iterations","option")
res2$ct_col = ct_col
res2$ct_col = factor(res2$ct_col, levels = ct_col)

res$iterations=res2$iterations


## Correlation
res$mean_corr=NA

for(i in 1:nrow(res)){
  qc = res$option[i]
  ct = res$cell_type[i]
  filenames = system(paste0("ls ./option_",qc,"/",ct,"_mx*"), intern = T)
  if(length(filenames)==0){next}
  data = read.table(filenames,header=T)
  
  corr = cor(data[,10:19])
  
  res$mean_corr[i] = mean( abs(as.numeric(corr[lower.tri(corr)])) )
  # res$median_corr[i] = median( abs(as.numeric(corr[lower.tri(corr)])) )
  
}


res$option=factor(res$option,levels=6:13)

pA <- ggplot(res, aes(x=option,y=mean_corr,group=cell_type)) +
  geom_line(linetype="dashed",alpha=0.1,col="grey") +
  geom_point(aes(col=cell_type,shape=cell_type),size=3) + 
  scale_y_continuous(
    name = "Mean pair-wise correlation"
  ) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=1:14) +
  scale_color_manual(values=ct_col) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.35,1,0.35,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "bottom") +
  xlab(" ") +
  scale_x_discrete(breaks=6:13,
                   labels=c(
                            "STD\npi1 excluded","STD\npi09 excluded","RINT\npi1 excluded","RINT\npi09 excluded",
                            "Log1p+STD\npi1 excluded","Log1p+STD\npi09 excluded",
                            "Log1p+STD\npi09 excluded\nHVG2000","RINT\npi09 excluded\nHVG2000") ) +
  guides(col=guide_legend(ncol=1),shape=guide_legend(ncol=1))
# facet_wrap(~cell_type,ncol=4,scales = "free")

# pA

pB <- ggplot(res, aes(x=option,y=iterations,group=cell_type)) +
  geom_line(linetype="dashed",alpha=0.1,col="grey") +
  geom_point(aes(col=cell_type,shape=cell_type),size=3) + 
  scale_y_continuous(
    name = "Iterations"
  ) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=1:14) +
  scale_color_manual(values=ct_col) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.35,1,0.35,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "right") +
  xlab(" ") +
  # scale_x_discrete(breaks=case,
  #                  labels=c("Not scaled","Not Standardization\npi1 excluded","Not Standardization\npi09 excluded","Log1p Standardization\npi1 excluded","Log1p Standardization\npi09 excluded","Standardization\npi1 excluded","Standardization\npi09 excluded")  ) +
  theme(axis.title.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank()) 
# facet_wrap(~cell_type,ncol=4,scales = "free")

# pB


pC<-ggplot(res, aes(x=option,y=time,group=cell_type)) +
  geom_line(linetype="dashed",alpha=0.1,col="grey") +
  geom_point(aes(col=cell_type,shape=cell_type),size=3) + 
  scale_y_continuous(
    name = "Time (mins)"
  ) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=1:14) +
  scale_color_manual(values=ct_col) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.35,1,0.35,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "bottom") +
  xlab("") +
  # scale_x_discrete(breaks=case,
  #                  labels=c("Not scaled","Not Standardization\npi1 excluded","Not Standardization\npi09 excluded","Log1p Standardization\npi1 excluded","Log1p Standardization\npi09 excluded","Standardization\npi1 excluded","Standardization\npi09 excluded")  ) +
  theme(axis.title.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
  guides(col=guide_legend(ncol=14),shape=guide_legend(ncol=1))
# facet_wrap(~cell_type,ncol=4,scales = "free")

# pC


# Combine plot
library(ggpubr)

t <- ggarrange(pA,pB,pC,ncol=1, nrow=3, common.legend = T,legend="right",labels=LETTERS[1:3],align ="v")
# t

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure5_time_iteration_mean_corr_v5.png",t,width=10,height=10,bg="white")


####