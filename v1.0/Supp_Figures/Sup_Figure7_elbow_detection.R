####
library(ggplot2)
library(ggpubr)
library(dplyr)

new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B","NK","NK_R","Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col=c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265","#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")

nr_PF = c()

for(k in 1:14){
alpha2=read.table(paste0("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/29APR2022/no_outliers/option_11/",new_names[k],"_variance_pf50.txt"),header=T)

x <- 1:nrow(alpha2)
y <- alpha2$Variance
x1 <- x[1]
y1 <- y[1]
x2 <- x[length(x)]
y2 <- y[length(y)]
x0 <- x
y0 <- y
distancesDenominator <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
distancesNumerator <- abs((x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1))
distances <- distancesNumerator/distancesDenominator
numOfPCsChosen <- which.max(distances)

nr_PF[k] = print(numOfPCsChosen)
}

nr_PF = as.data.frame(cbind(ct = new_names, elbow = nr_PF))
nr_PF$elbow = as.numeric(nr_PF$elbow)

####
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

# Calculate the increase % of PFs fitted
old$percent = NA
for(i in 1:nrow(old)){
  if(old$PF[i]==0){next}
  old$percent[i] = (old$eGene[i] - old$eGene[i-1]) / old$eGene[i-1]
}
# Fit smooth curve
nr = c()
for(i in 1:14){
  
  tmp = old[old$ct == new_names[i] & old$Group=="All_genes",]
  tmp = tmp[-1,]
  tmp <- tmp %>% 
    mutate(smooth_y2 = predict(loess(percent ~ PF, span = 0.15, data = tmp)))
  
  for(k in 1:50){
    if(tmp$smooth_y2[k]<0){break}
  }
  nr[i] = print(k - 1)
  
}

nr_PF$local_greedy = nr

# Merge results
old = old[old$Group == "All_genes",]
old = merge(old, nr_PF, by = "ct", sort = F)

#### eQTL sensitivity plot with optimal K chosen
p7 <- ggplot(old, aes(x = PF, y = eGene)) +
  geom_point(aes(shape = Group, col = ct)) + 
  geom_vline(aes(xintercept = elbow)) +
  geom_vline(aes(xintercept = local_greedy), linetype = "dashed") +
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
  xlab("Nr of PEER factors") +
  theme(axis.title.x = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(vjust = -1),axis.title.y = element_text(vjust = 1)) +
  facet_wrap(~ct,ncol=4,scales = "free") +
  guides(col=guide_legend(ncol=3))

p7

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/Sup_Figure7_elbow_detection.png", p7, width = 10.5, height = 8, bg = "white")

# Power gain % comparison
pw = data.frame(ct = nr_PF$ct, power_gain = NA, power_gain_elbow = NA)
for(i in 1:14){
  x1=old[old$ct==nr_PF$ct[i] & old$PF == 0, "eGene"]
  x2=old[old$ct==nr_PF$ct[i] & old$PF == nr_PF$local_greedy[i], "eGene"]
  pw$power_gain[i] = (x2 - x1) / x1
  
  x3=old[old$ct==nr_PF$ct[i] & old$PF == nr_PF$elbow[i], "eGene"]
  pw$power_gain_elbow[i] = (x3 - x1) / x1
  }


t1 = read.table("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Table/Table_S1_eGene_power_gain.txt",header=T)
pw = pw[match(t1$ct,pw$ct),]

