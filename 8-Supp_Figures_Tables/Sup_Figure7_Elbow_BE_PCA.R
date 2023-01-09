####
library(ggplot2)

new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
               "Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col=c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987",
         "#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")

res = c()

for(k in 1:13){
  a = c()
for(i in 1:14){
  b = read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/PCA/Elbow_BE/option_",k,"/",new_names[i],"_Elbow_BE_pca50.txt"),header=T)
  b$ct = new_names[i]
  b$ct_col = ct_col[i]
  b$option = k
  a=rbind(a,b)
}

res = rbind(res,a)

}

# 
res$ct = factor(res$ct, levels = new_names[1:14])

res1=res[,-2]
res2=res[,-1]
colnames(res1)[1] = colnames(res2)[1] = "Nr_PCs"
res1$Method = "Elbow"
res2$Method = "BE"
res = rbind(res1,res2)

p12 <- ggplot(res, aes(x = option, y = Nr_PCs)) +
  # geom_smooth(method = "loess", se = F, fullrange = FALSE, level = 0.95, col = "black", aes(linetype = Method)) +
  geom_point(aes(shape = Method, col = ct)) + 
  # geom_vline(aes(xintercept = nr_PC)) +
  # geom_vline(aes(xintercept = nr_PF), linetype = "dashed") +
  scale_y_continuous(
    name = "Nr of PCs"
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
  xlab("QC option #") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  scale_x_continuous(breaks=c(1:13),
                   labels=c(1:13)) + 
  facet_wrap(~ct,ncol=4,scales = "free") +
  guides(col=guide_legend(ncol=3))

# p11 <- ggarrange(p11a, p11b, ncol = 1, nrow = 2, common.legend = F, labels=LETTERS[1:2], align = "v")

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure7_Elbow_BE_pca50.png", p12, width = 9.5, height = 7, bg = "white")



