#### Make Figure 2  ####
library(ggplot2)
library(ggpubr)

# new_names <- c("B_IN", "B_MEM", "CD4_NC", "CD4_ET", "CD4_SOX4", "CD8_NC", "CD8_ET",
#                "CD8_S100B", "DC", "Mono_C", "Mono_NC", "NK_R", "NK", "Plasma",
#                "Erythrocytes","Platelets")

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
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_raw_PF1_50.txt"),header=T,sep="\t")
  # b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/new/17APR2022/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  b$ct_col=ct_col[i]
  a=rbind(a,b)
}

new=a

old$Group="QC option #11"
new$Group="No QC"

old=rbind(old,new)

old$ct = factor(old$ct, levels = new_names[1:14])
old$Group = factor(old$Group, levels = c("QC option #11","No QC"))

p1<-ggplot(old, aes(x=PF,y=eQTL, group=Group)) +
  #geom_abline(slope=1,intercept=0,lty=3) +
  geom_smooth(method="loess", col="black", se=T, fullrange=FALSE, level=0.95, aes(linetype = Group)) +
  geom_point(aes(shape=Group,col=ct)) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Nr of eQTLs"
  ) +
  scale_color_manual(values=ct_col) +
  scale_shape_manual(values=c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.box = "horizontal",
        legend.position = c(1, 0),
        legend.justification = c(1, 0)) +
  xlab("Nr of PEER factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct,ncol=4,scales = "free") +
  guides(col=guide_legend(ncol=3))

# p1
p2<-ggplot(old, aes(x=PF,y=eGene, group=Group)) +
  #geom_abline(slope=1,intercept=0,lty=3) +
  geom_smooth(method="loess", col="black", se=T, fullrange=FALSE, level=0.95, aes(linetype = Group)) +
  geom_point(aes(shape=Group,col=ct)) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Nr of eGenes"
  ) +
  scale_color_manual(values=ct_col) +
  scale_shape_manual(values=c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.box = "horizontal",
        legend.position = c(1, 0),
        legend.justification = c(1, 0)) +
  xlab("Nr of PEER factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct,ncol=4,scales = "free") +
  guides(col=guide_legend(ncol=3))


t <- ggarrange(p1,p2,ncol=1, nrow=2, labels=LETTERS[1:2])

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure5_nr_eGene_no_QC_vs_QC_option_11_v2.png",t,width=10.5,height=15,bg="white")


## % increase of eGene detection power
library(dplyr)

inc = old %>% group_by(ct,Group) %>% filter(eGene == max(eGene))
inc = as.data.frame(inc)

inc = inc[!duplicated(inc$eGene),]

summary((inc[1:14,"eGene"] - inc[15:28,"eGene"]) / (inc[15:28,"eGene"]))
summary((inc[1:14,"eQTL"] - inc[15:28,"eQTL"]) / (inc[15:28,"eQTL"]))

