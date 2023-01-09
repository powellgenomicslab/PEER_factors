####
library(corrplot)

pdf(file = "~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure6_PF_cor_Drew_iPSC.pdf",width=10,height = 8)
par(mfrow=c(2,2),mar=c(4,5,5,1),xpd=T)

for(i in 1:4){
  f1=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/cluster",i,"_iPSC_peer_factors.txt"),header=T)
  colnames(f1)[-1]=paste0("PF",1:10)
  M = cor(f1[,-1])
  corrplot.mixed(M,title=paste("iPSC Cluster",i),mar=c(0,0,4.5,0),xpd=T) # colorful number
  # corrplot(M, title=paste("Cluster",i),method = 'circle', type = 'lower', insig='blank',
  #          addCoef.col ='black', number.cex = 0.8, diag=FALSE,mar=c(0,0,1,0))
}
dev.off()

# Fibroblast
pdf(file = "~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure6_PF_cor_Drew_fibroblast.pdf",width=10,height = 12)
par(mfrow=c(3,2),mar=c(4,5,5,1),xpd=T)

for(i in 1:6){
  f1=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/cluster",i,"_peer_factors.txt"),header=T)
  colnames(f1)=paste0("PF",1:10)
  M = cor(f1)
  corrplot.mixed(M,title=paste("Fibroblast Cluster",i),mar=c(0,0,4.5,0),xpd=T) # colorful number
  # corrplot(M, title=paste("Cluster",i),method = 'circle', type = 'lower', insig='blank',
  #          addCoef.col ='black', number.cex = 0.8, diag=FALSE,mar=c(0,0,1,0))
}
dev.off()

#### Down-sampling fibroblast clusters ####

pw_cor=c()
for(i in 1:6){
tmp=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/down31/stderr_PEER_fibroblast_cluster",i),header=F)
colnames(tmp)="cor"
tmp$cluster=paste0("Cluster",i)
tmp$N=31

tmp2=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/down40/stderr_PEER_fibroblast_cluster",i),header=F)
colnames(tmp2)="cor"
tmp2$cluster=paste0("Cluster",i)
tmp2$N=40

tmp3=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/down50/stderr_PEER_fibroblast_cluster",i),header=F)
colnames(tmp3)="cor"
tmp3$cluster=paste0("Cluster",i)
tmp3$N=50


pw_cor=rbind(pw_cor,tmp, tmp2, tmp3)
}

p_cor <- ggplot(pw_cor, aes(x=N, y=cor, fill=N, group=N)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + geom_jitter(width = 0.2, size = 1, pch=21) +
  ylab("Mean of pair-wise correlation") +
  xlab("Sample size (N)") +
  facet_wrap(~cluster,nrow = 3, scales = "free")

p_cor

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure7_downsampling_v2.png",p_cor,width=8,height=8)


#### Each panel corrplot ####
i=1
f1=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/cluster",i,"_iPSC_peer_factors.txt"),header=T)
pdf(file="~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/pair_covar_corplot_iPSC_cluster1.pdf",width = 9,height = 7)
pairs.panels(f1[,-1], cex.cor=1)
dev.off()

f1=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/cluster",i,"_peer_factors.txt"),header=T)
pdf(file="~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/Drew/pair_covar_corplot_fibroblast_cluster1.pdf",width = 9,height = 7)
pairs.panels(f1[,-1], cex.cor=1)
dev.off()