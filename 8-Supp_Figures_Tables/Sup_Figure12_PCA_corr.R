####
library(ggplot2)
library(ggpubr)

setwd("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/PCA/")
new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R",
               "Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col=c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987",
         "#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")

res = c()

for(i in 1:13){
  tmp = data.frame(ct = new_names, option = i, mean_corr = NA, median_corr = NA, max_corr = NA)
  
  for(k in 1:14){
    
    data = read.table(paste0("./option_", i, "/", new_names[k], "_mx_pca50.txt"), header = T)
    corr = cor(data[,10:59])
    
    tmp$mean_corr[k] = mean( abs(as.numeric(corr[lower.tri(corr)])) )
    tmp$median_corr[k] = median( abs(as.numeric(corr[lower.tri(corr)])) )
    tmp$max_corr[k] = max( abs(as.numeric(corr[lower.tri(corr)])) )
  }
  
  res = rbind(res, tmp)
  print(i)
}

res$option=factor(res$option,levels=1:13)


# plot
p10a <- ggplot(res, aes(x = option, y = mean_corr, group = ct)) +
  geom_line(linetype="dashed",alpha=0.3,col="grey") +
  geom_point(aes(col = ct,shape = ct),size=3) + 
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
  guides(col = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))

p10b <- ggplot(res, aes(x = option, y = max_corr, group = ct)) +
  geom_line(linetype="dashed",alpha=0.3,col="grey") +
  geom_point(aes(col = ct,shape = ct),size=3) + 
  labs(y = "Max pair-wise correlation") +
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
  guides(col = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))


p10 <- ggarrange(p10a, p10b, ncol = 1, nrow = 2, common.legend = T, legend = "bottom", labels=LETTERS[1:2], align = "v")

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure12_mean_max_corr_PCA.png", p10, width = 13, height = 6, bg = "white")


