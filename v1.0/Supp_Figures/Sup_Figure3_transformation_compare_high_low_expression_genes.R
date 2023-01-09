#### Why log(x+1) and RINT is not good enough?
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(RNOmni)
library(cowplot)

setwd("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/presentation/")

# Real data
mean_mx=fread("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/mean_var_matrix/no_outliers/CD4_NC_cells_mean_mx.txt",header=T,check.names=FALSE)
mean_mx=as.data.frame(mean_mx)
gene=fread("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/mean_var_matrix/no_outliers/CD4_NC_cells_gene_list.txt",header=T,check.names=FALSE)
gene=as.data.frame(gene)
mean_mx=mean_mx[,which(!is.na(mean_mx[1,]))]

m=as.numeric(rowMeans(mean_mx))
sk=apply(mean_mx,1,function(x)skewness(x))

# High expression, FTH1
g_symbol="FTH1"
x1=as.numeric(mean_mx[which(gene$gene==g_symbol),])

df1 <- data.frame(x = x1, y = log1p(x1))
p1 <- ggplot(df, aes(x, y)) + 
  geom_point() + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, vjust=-1, face = "italic")) +
  ggtitle(g_symbol) +
  geom_abline(slope=1,intercept=0,lty=3,col="grey") + 
  expand_limits(x = 0, y = 0) +
  xlab("Expression") +
  ylab("Log(x+1) transformed")

p1 <- ggExtra::ggMarginal(p1, type = "histogram")
# p1

# RINT
x2=as.numeric(mean_mx[which(gene$gene==g_symbol),])

df2 <- data.frame(x = x2, y = RankNorm(x2))
p2 <- ggplot(df2, aes(x, y)) + 
  geom_point() + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, vjust=-1)) +
  ggtitle(" ") +
  # geom_abline(slope=1,intercept=0,lty=3,col="grey") + 
  expand_limits(x = 0, y = 0) +
  xlab("Expression") +
  ylab("RINT transformed")

p2 <- ggExtra::ggMarginal(p2, type = "histogram")
# p2


# ggsave(paste0("real_data_CD4_NC_",g_symbol,"_log1p.pdf"),plot=p4,width=4,height=3.5)

# Moderate expression, CDC42

# Extremely low expression, PMCH
g_symbol="PMCH"
x3=as.numeric(mean_mx[which(gene$gene==g_symbol),])

df3 <- data.frame(x = x3, y = log1p(x3))
p3 <- ggplot(df3, aes(x, y)) + 
  geom_point() + 
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0, vjust=-1, face = "italic")) +
  ggtitle(g_symbol) +
  geom_abline(slope=1,intercept=0,lty=3,col="grey") + 
  expand_limits(x = 0, y = 0) +
  xlab("Expression") +
  ylab("Log(x+1) transformed")

p3 <- ggExtra::ggMarginal(p3, type = "histogram")
# p3

# ggsave(paste0("real_data_CD4_NC_",g_symbol,"_log1p.pdf"),plot=p4,width=4,height=3.5)
# RINT
x4=as.numeric(mean_mx[which(gene$gene==g_symbol),])

df4 <- data.frame(x = x4, y = RankNorm(x4))
p4 <- ggplot(df4, aes(x, y)) + 
  geom_point() + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0, vjust=-1)) +
  ggtitle(" ") +
  # geom_abline(slope=1,intercept=0,lty=3,col="grey") + 
  expand_limits(x = 0, y = 0) +
  xlab("Expression") +
  ylab("RINT transformed")

p4 <- ggExtra::ggMarginal(p4, type = "histogram")
# p4

# t1 <- ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = F,legend="none",labels=LETTERS[1:4])


t1 <- ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(p3, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_plot(p4, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5))

ggsave(filename = "~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure3_transformation_compare.pdf",t1,width=10,height=9)

####


####