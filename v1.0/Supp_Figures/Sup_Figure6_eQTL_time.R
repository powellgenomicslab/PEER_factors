####
library(readr)

setwd("~/Documents/Study/Post-doc/project/vQTL/results/MatrixQTL/PF_test/time/")

new_names <- c("CD4_NC","CD4_ET","CD4_SOX4","CD8_ET","CD8_NC","CD8_S100B",
               "NK","NK_R","Plasma","B_MEM","B_IN","Mono_C","Mono_NC","DC")

ct_col = c("#882E72","#B178A6","#D6C1DE","#1965B0","#5289C7","#7BAFDE","#4EB265",
         "#90C987","#CAE0AB","#F7EE55","#F6C141","#F1932D","#E8601C","#DC050C")

res = c()

for(k in 1:14){

tmp = c()

for(i in 0:50){

a = read.table(paste0("./time_", new_names[k], "_PF", i,".txt"), header = F)
a = a[,c(1,4)]
a$V1 = parse_number(a$V1)
colnames(a) = c("CHR", "Time")
a$PF = i

tmp = rbind(tmp, a)
}

tmp$ct = new_names[k]
tmp$ct_col = ct_col[k]
res = rbind(res, tmp)
print(new_names[k])

}

res$Time = res$Time / 60
res$ct = factor(res$ct, levels = new_names[1:14])

p6 <- ggplot(res, aes(y = Time)) +
  #geom_abline(slope=1,intercept=0,lty=3) +
  geom_boxplot(aes(x = as.factor(PF), col = ct), outlier.size = 0.05, width = 0.4, lwd = 0.2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Time (mins)"
  ) +
  scale_color_manual(values = ct_col) +
  # scale_shape_manual(values=c(1,2)) +
  theme_minimal() +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  theme(axis.title.y = element_text(color = "black", size=13)) +
  theme(legend.title= element_blank()) +
  theme(legend.box = "horizontal",
        legend.position = c(0.95, 0),
        legend.justification = c(1, 0)) +
  xlab("Nr of PEER factors") +
  # theme(axis.title.x = element_text()) +
  theme(axis.title.x = element_text(hjust = 0.5, vjust = -1, angle = 90, size = rel(0.2)), axis.title.y = element_text(vjust = 1)) +
  scale_x_discrete(breaks = seq(0, 50, 5),
                   labels = seq(0, 50, 5)) +
  facet_wrap(~ct, ncol = 4) +
  guides(col = guide_legend(ncol = 3))

# p6

ggsave("~/Documents/Study/Post-doc/project/vQTL/manuscipt/PF_commentary/Sup_Figure/FigureS7_eQTL_time.png", 
       p6, width = 9.5, height = 7, bg = "white")

