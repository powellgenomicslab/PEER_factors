####
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(ggExtra)
library(GGally)
library(psych)
library(viridis)
library(cowplot)
library(ggpubr)

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

pA = list()

for(k in 1:14){

## Figure 1A 
old=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_1/",new_names[k],"_mx_pf50.txt"),header=T,check.names = F)

# id=old$id
# old=old[,-1]
# old=t(old)
# old=cbind(rownames(old),old)
# colnames(old)=c("id",id)
# old=as.data.frame(old)
# rownames(old)=1:nrow(old)

old[,-1] = as.data.frame(sapply(old[,-1], as.numeric))
colnames(old)=c("ID","Sex",paste0("PC",1:6),"Age",paste0("PF",1:50))

par(mar=c(0,0,0,0))
pairs.panels(old[,10:19], 
             hist.col = viridis::viridis(10)[4],
             smoother = F,
             smooth = F,
             stars=T,
             ellipses = F, # show correlation ellipses
             cex = 0.1,
             cex.labels = 0.9,
             cex.cor = 2, xaxt = 'n', yaxt = 'n', alpha=0.02)

pA[[k]] <- recordPlot()

dev.off(dev.list()["RStudioGD"])

#### Figure 1B: Variance Elbow plot of PFs ####

# alpha1=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/raw_SCT/pf50/",new_names2[k],"_variance_no_scale_pi1_all_genes_PF50.txt"),header=T)
alpha1=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_1/",new_names[k],"_variance_pf50.txt"),header=T)
alpha1$option="Raw"
alpha2=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_11/",new_names[k],"_variance_pf50.txt"),header=T)
alpha2$option="QC"

alpha=rbind(alpha1[1:10,],alpha2[1:10,])

# ggplot
pB1<-ggplot(alpha[alpha$option=="Raw",], aes(x=1:10,y=Variance)) +
  geom_line(linetype="dashed",alpha=0.9,col=viridis::viridis(10)[4]) +
  geom_point(aes(size=3),col=viridis::viridis(10)[4]) + 
  theme_minimal() +
  # theme(plot.margin = unit(c(0.35,1,0.35,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "none") +
  xlab("Factors") +
  ylab("Weights") +
  scale_x_discrete(name="Factors",limits=factor(1:10),labels=1:10, breaks=1:10, guide = guide_axis(check.overlap = TRUE))

pB2<-ggplot(alpha[alpha$option=="QC",], aes(x=1:10,y=Variance)) +
  geom_line(linetype="dashed",alpha=0.9,col=viridis::viridis(10)[1]) +
  geom_point(aes(size=3),col=viridis::viridis(10)[1]) + 
  theme_minimal() +
  # theme(plot.margin = unit(c(0.35,1,0.35,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.position = "none") +
  xlab("Factors") +
  ylab("Weights") +
  scale_x_discrete(name="Factors",limits=factor(1:10),labels=1:10, breaks=1:10, guide = guide_axis(check.overlap = TRUE))

pB <- ggarrange(pB1,pB2,ncol=1, nrow=2, common.legend = T,legend="none",align = "v")

#### Figure 1C: Per gene mean against Fano factor ####
# tryCatch(mean_mx=fread(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/mean_var_matrix/no_outliers/",new_names[k],"_cells_mean_mx.txt"),header=T,check.names=FALSE), finally=freadCleanup)
mean_mx=fread(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/mean_var_matrix/no_outliers/",new_names[k],"_cells_mean_mx.txt"),header=T,check.names=FALSE)
mean_mx=as.data.frame(mean_mx)

mean_mx=mean_mx[,which(!is.na(mean_mx[1,]))]

m=rowMeans(mean_mx,na.rm=T)
v=apply(mean_mx,1,function(x)var(x,na.rm=T))

data=as.data.frame(cbind(m,v))
data$ff=data$v/data$m
data$pi=rowSums(mean_mx==0,na.rm=T)/sum(!is.na(mean_mx[1,]))

colorsP <- brewer.pal(10, "Spectral")
labP <- c("0;0.1","0.1;0.2","0.2;0.3","0.3;0.4","0.4;0.5","0.5;0.6","0.6;0.7","0.7;0.8","0.8;0.9","0.9;1.0")
vecP <- seq(0,1,0.1)
data$col_group=cut(data$pi,breaks=vecP,labels=labP,include.lowest = T, right = T)

# Skewness
# require(moments)
# data$skewness=apply(mean_mx,1,function(x)skewness(x))

## Mean vs Fano factor
pc<-ggplot(data, aes(x=log10(m),y=log10(ff))) +
  # geom_tile(aes(fill=col_group),size=1) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  geom_point(aes(col = pi), size = 1, alpha = 0.5) +
  theme(legend.position="bottom",
        legend.text = element_text(angle = 0)) +
  #scale_color_brewer(palette = "Spectral") +
  labs(x=expression(paste(Log[10],"(mean)",sep="")),
       y=expression(paste(Log[10],"(Fano factor)",sep=""))) +
  scale_colour_binned(name = expression(Proportion~of~zeros~(pi[0])~"  "),
                      type= "viridis",
                      breaks = vecP,
                      limits = c(0, 1),
                      guide = guide_coloursteps(even.steps = FALSE,
                                                show.limits = F,
                                                barwidth=15)) 

pC <- ggExtra::ggMarginal(pc, type = "histogram",
                          fill = viridis::viridis(10)[4])

#### Figure 1D: compare PFs between all genes and top 2000 HVGs ####
a=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_11/",new_names[k],"_mx_pf50.txt"),header=T)
b=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_12/",new_names[k],"_mx_pf50.txt"),header=T)

res=c()
for(i in 1:10){
  tmp=a[,9+i]
  tmp=as.data.frame(tmp)
  tmp$PF=paste0("PF",i)
  res=rbind(res,tmp)
}

res2=c()
for(i in 1:10){
  tmp=b[,9+i]
  tmp=as.data.frame(tmp)
  tmp$PF=paste0("PF",i)
  res2=rbind(res2,tmp)
}

colnames(res)[1]="All_genes"
res$HVG2000=res2$tmp

res$PF = factor(res$PF, levels = paste0("PF",1:10))

pD<-ggplot(res, aes(x=All_genes,y=HVG2000)) +
  # geom_tile(aes(fill=col_group),size=1) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text=element_text(size=6)) +
  geom_point(aes(col=PF),size=1,alpha=0.3) +
  theme(legend.position="none",
        legend.text = element_text(angle = 0)) +
  facet_wrap(~PF,nrow = 3,scales = "free") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  scale_color_brewer(palette = "Paired")
  # scale_color_viridis(discrete=TRUE)

## Combine the four plots
# pdf(file = paste0("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure1/Panel_",new_names[k],".pdf"),width=13,height=12,bg="white")

# png(filename = paste0("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure1/Panel_",new_names[k],".png"), width=13*120, height=12*120, bg = "white",type = "cairo")

pall <- ggdraw() +
  draw_plot(pA[[k]], x = 0, y = .45, width = .7, height = .55) +
  draw_plot(pB, x = .69, y = .45, width = .3, height = .55) +
  draw_plot(pC, x = 0, y = 0, width = 0.45, height = 0.45) +
  draw_plot(pD, x = 0.48, y = 0, width = 0.5, height = 0.45) +
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0, 0.69, 0, 0.45), y = c(1, 1, 0.45, 0.45)) +
  draw_label(new_names3[k], fontface='bold', x = 0.35, y = 0.98)

# ggsave(filename = paste0("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure1/Panel_",new_names[k],"_v2.pdf"),plot=pall,width=13,height=12,bg="white")
ggsave(filename = paste0("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure1/Panel_",new_names[k],"_v2.png"),plot=pall,width=13,height=12,bg="white")

# dev.off()

print(new_names[k])

rm(list=setdiff(ls(), c("new_names","new_names2","new_names3","pA")))

}
