####
library(ggplot2)

new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
               "Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col=c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987",
         "#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")

a=c()

for(i in 1:14){
  #b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/PCA/all_genes/",new_names[i],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_all_genes_PC1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  b$ct_col=ct_col[i]
  a=rbind(a,b)
}

old=a

a=c()

for(i in 1:14){
  b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_all_genes_PF1_50.txt"),header=T,sep="\t")
  # b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/new/17APR2022/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
  b$ct=new_names[i]
  b$ct_col=ct_col[i]
  a=rbind(a,b)
}

new=a
# 
# old = old[,-c(4)]
# new = new[,-c(4,5)]

# 
old$Group="PCA_all_genes"
new$Group="PEER_all_genes"

old=rbind(old,new)

old$ct = factor(old$ct, levels = new_names[1:14])

p11a <- ggplot(old1, aes(x = PF, y = eGene)) +
  geom_point(aes(shape = Group, col = ct)) + 
  # geom_vline(aes(xintercept = nr_PC)) +
  # geom_vline(aes(xintercept = nr_PF), linetype = "dashed") +
  geom_smooth(method = "loess", se = F, fullrange = FALSE, level = 0.95, col = "black", aes(linetype = Group)) +
  scale_y_continuous(
    # Features of the first axis
    name = "Nr of eGene"
  ) +
  scale_color_manual(values = ct_col) +
  scale_shape_manual(values = c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.box = "horizontal",
        legend.position = c(0.95, 0),
        legend.justification = c(1, 0)) +
  xlab("Nr of latent factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct,ncol=4,scales = "free") +
  guides(col=guide_legend(ncol=3))


## HVG2000
a=c()

for(i in 1:14){
  #b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/",new_names[i],"_number_eQTL_signals_mx_RINT_pi09_PF1_50.txt"),header=T,sep="\t")
    b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/PCA/HVG2000/",new_names[i],"_number_eQTL_signals_mx_log1p_std_pi09_excluded_HVG2000_PC1_50.txt"),header=T,sep="\t")
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
# 
old = old[,-c(4)]
new = new[,-c(4,5)]

# 
old$Group="PCA_HVG2000"
new$Group="PEER_HVG2000"

old=rbind(old,new)

old$ct = factor(old$ct, levels = new_names[1:14])

p11b <- ggplot(old, aes(x = PF, y = eGene)) +
  geom_point(aes(shape = Group, col = ct)) + 
  # geom_vline(aes(xintercept = nr_PC)) +
  # geom_vline(aes(xintercept = nr_PF), linetype = "dashed") +
  geom_smooth(method = "loess", se = F, fullrange = FALSE, level = 0.95, col = "black", aes(linetype = Group,size=0.5)) +
  scale_y_continuous(
    # Features of the first axis
    name = "Nr of eGene"
  ) +
  scale_color_manual(values = ct_col) +
  scale_shape_manual(values = c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.box = "horizontal",
        legend.position = c(0.95, 0),
        legend.justification = c(1, 0)) +
  xlab("Nr of latent factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct,ncol=4,scales = "free") +
  guides(col=guide_legend(ncol=3))

p11 <- ggarrange(p11a, p11b, ncol = 1, nrow = 2, common.legend = F, labels=LETTERS[1:2], align = "v")

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure13_PCA_eQTL_sensitivity_test.png", p11, width = 13, height = 16, bg = "white")



