# compare GC/AT content in 3'UTR between aOutlier genes and eOutlier genes

# load data
dat.m <- readRDS("GC_AT_content.aOutlier_vs_eOutlier.RDS")
library(cowplot)
library (ggplot2)
library(ggsci)
library(ggpubr)

# plot
p1 <- ggplot(dat.m,aes(x=Group,y=ratio,fill=Group)) + geom_boxplot(width=.5) + theme_pubr() + 
  scale_fill_manual(breaks = c("eoutlier","aoutlier"),values = c("#4DBBD5FF","#E64B35FF"))

p2<- ggplot(dat.m,aes(x=Group,y=at_ratio,fill=Group)) + geom_boxplot(width=.5) + theme_pubr() + 
  scale_fill_manual(breaks = c("eoutlier","aoutlier"),values = c("#4DBBD5FF","#E64B35FF")) +
  ylab("ratio")

# test the difference
wilcox.test(aoutlier_dat_ratio$ratio,eoutlier_dat_ratio$ratio,alternative = "less")
wilcox.test(aoutlier_dat_ratio$at_ratio,eoutlier_dat_ratio$at_ratio,alternative = "greater")
