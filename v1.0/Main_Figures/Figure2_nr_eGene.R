#### Make Figure 2  ####
library(ggplot2)
library(tidyr)
library(readr)

new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R","Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

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

new = a

#
old = old[,-c(4)]
new = new[,-c(4,5)]

old$Group = "All_genes"
new$Group = "HVG2000"

old = rbind(old,new)

old$ct = factor(old$ct, levels = new_names[1:14])

p1 <- ggplot(old, aes(x = PF, y = eGene)) +
  #geom_abline(slope=1,intercept=0,lty=3) +
  geom_point(aes(shape = Group, col = ct)) + 
  # geom_smooth(method="loess", se=T, fullrange=FALSE, level=0.95, aes(linetype = Group)) +
  scale_y_continuous(
    # Features of the first axis
    name = "Nr of eGene"
  ) +
  scale_color_manual(values = ct_col) +
  scale_shape_manual(values = c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title = element_blank()) +
  theme(legend.box = "horizontal",
        legend.position = c(0.95, 0),
        legend.justification = c(1, 0)) +
  xlab("Nr of PEER factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1), axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct, ncol = 4, scales = "free") +
  guides(col = guide_legend(ncol = 3))

#### Panel for time and iterations
setwd("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/")

res=c()

for(i in 1:13){
  
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

res$cell_type = factor(res$cell_type,levels = new_names)

res$ct_col = ct_col
res$ct_col = factor(res$ct_col, levels = ct_col)

pB<-ggplot(res, aes(x=option,y=time,group=cell_type)) +
  geom_line(linetype="dashed",alpha=0.3,col="grey") +
  geom_point(aes(col=cell_type,shape=cell_type),size=3) + 
  scale_y_continuous(
    name = "Time (mins)"
  ) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=1:14) +
  scale_color_manual(values=ct_col) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.2,1,0.2,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "bottom") +
  xlab("") +
  # scale_x_discrete(breaks=case,
  #                  labels=c("Not scaled","Not Standardization\npi1 excluded","Not Standardization\npi09 excluded","Log1p Standardization\npi1 excluded","Log1p Standardization\npi09 excluded","Standardization\npi1 excluded","Standardization\npi09 excluded")  ) +
  theme(axis.title.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
  guides(col=guide_legend(ncol=14),shape=guide_legend(nrow=1))
# facet_wrap(~cell_type,ncol=4,scales = "free")

#### Iterations
res2=c()

for(i in 1:13){
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

pB2 <- ggplot(res, aes(x=option,y=iterations,group=cell_type)) +
  geom_line(linetype="dashed",alpha=0.3,col="grey") +
  geom_point(aes(col=cell_type,shape=cell_type),size=3) + 
  scale_y_continuous(
    name = "Iterations"
  ) +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=1:14) +
  scale_color_manual(values=ct_col) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.2,1,0.2,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "right") +
  xlab(" ") +
  # scale_x_discrete(breaks=case,
  #                  labels=c("Not scaled","Not Standardization\npi1 excluded","Not Standardization\npi09 excluded","Log1p Standardization\npi1 excluded","Log1p Standardization\npi09 excluded","Standardization\npi1 excluded","Standardization\npi09 excluded")  ) +
  theme(axis.title.x = element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank()) 
# facet_wrap(~cell_type,ncol=4,scales = "free")

# p2

## Correlation
res$mean_corr=NA
res$median_corr=NA

for(i in 1:nrow(res)){
  qc = res$option[i]
  ct = res$cell_type[i]
  filenames = system(paste0("ls ./option_",qc,"/",ct,"_mx*"), intern = T)
  if(length(filenames)==0){next}
  data = read.table(filenames,header=T)
  
  corr = cor(data[,10:59])
  
  res$mean_corr[i] = mean( abs(as.numeric(corr[lower.tri(corr)])) )
  res$median_corr[i] = median( abs(as.numeric(corr[lower.tri(corr)])) )
  
}

res$option=factor(res$option,levels=1:13)

pA <- ggplot(res, aes(x=option,y=mean_corr,group=cell_type)) +
  geom_line(linetype="dashed",alpha=0.3,col="grey") +
  geom_point(aes(col=cell_type,shape=cell_type),size=3) + 
  labs(y = "Mean pair-wise correlation") +
  # scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=1:14) +
  scale_color_manual(values=ct_col) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.5,1,0.2,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13, vjust = 0.1)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "bottom") +
  xlab(" ") +
  scale_x_discrete(breaks=1:13,
                   labels=c("None","pi1 excluded","pi09 excluded","Log1p\npi1 excluded","Log1p\npi09 excluded",
                            "STD\npi1 excluded","STD\npi09 excluded","RINT\npi1 excluded","RINT\npi09 excluded",
                            "Log1p+STD\npi1 excluded","Log1p+STD\npi09 excluded",
                            "Log1p+STD\npi09 excluded\nHVG2000","RINT\npi09 excluded\nHVG2000") ) +
  guides(col=guide_legend(nrow=1),shape=guide_legend(nrow=1))


# Combine plot
library(ggpubr)

t1 <- ggarrange(pA,pB,ncol=1, nrow=2, common.legend = T,legend="bottom",labels=LETTERS[1:2],align = "v")

t2 <- ggarrange(t1,p1,ncol=1, nrow=2, common.legend = F,labels=LETTERS[c(1,3)],heights=c(2.2,3))

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Figure/Figure2_mean_corr_time_eQTL.pdf",t2,width=13,height=14,bg="white")
ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Figure/Figure2_mean_corr_time_eQTL.png",t2,width=13,height=14,bg="white")


####

#### Find turning point ####
nr_peer=data.frame(ct=new_names[1:14],nr_peer=NA)
tmp=old[old$Group=="All_genes",]

for(i in 1:14){
  a=tmp[tmp$ct==new_names[i],]
  nr_peer$nr_peer[i]=which.max(a$eGene)
  print(which.max(a$eGene))
  
}

nr_peer[nr_peer$ct %in% c("CD4_NC","CD8_NC","CD8_ET","NK"),2]=20
tmp=tmp[tmp$ct=="CD8_S100B",]
tmp=tmp[tmp$PF<20,]
nr_peer[nr_peer$ct=="CD8_S100B",2]=which.max(tmp$eGene)

write.table(nr_peer,"~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/new/Nr_PEER_factors_fit_for_each_celltype.txt",row.names = F,col.names = T,quote = F)

#### New settings with local FDR ####
####
library(ggplot2)

new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
               "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
               "Erythrocytes","Platelets")

a=c()

for(i in 1:14){
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/new/17APR2022/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  a=rbind(a,b)
}

old=a

a=c()

for(i in 1:14){
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/new/17APR2022/HVG2000/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  a=rbind(a,b)
}

new=a

old$Group="All_genes"
new$Group="HVG2000"

old=rbind(old,new)

old$ct = factor(old$ct, levels = new_names[1:14])


p1<-ggplot(old, aes(x=PF)) +
  #geom_abline(slope=1,intercept=0,lty=3) +
  geom_point(aes(y=eGene,shape=Group,col=Group)) + 
  # col=c("#E69F00","#56B4E9")
  scale_y_continuous(
    # Features of the first axis
    name = "Nr of eGene"
  ) +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  scale_shape_manual(values=c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = c(0.81, 0.01),
        legend.justification = c(1, 0)) +
  xlab("Nr of PEER factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct,ncol=4)

p1

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Figure/Figure2_nr_eGene_localFDR_19APR2022.pdf",p1,width=9.5,height=7)








#### END