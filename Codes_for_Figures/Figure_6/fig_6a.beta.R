setwd("/Users/xudong.zou/SynologyDrive/Projects/Manuscripts/xxx/Figures/Figure6/")
rm(list=ls())


library(data.table)
library(dplyr)
library(magrittr)



dat <- fread("./data/Combined_matched_RVs.outlier_coloc.txt",header=F,sep="\t")
names(dat) <- c("chrom","pos1","posterior","MAF","isOutlier","isColoc","Beta_PCT")
dat %<>% filter(isOutlier!=1)
dat$Group <- "nonOutlier"
dat$Group[dat$isOutlier==2] <- "Outlier"
dat$Group[dat$isOutlier==2 & dat$isColoc==1] <- "Coloc-Outlier"

library(ggplot2)
library(ggpubr)
ggplot(dat,aes(x=Group,y=Beta_PCT)) + geom_boxplot(width=.5) + theme_pubr()
ggplot(dat,aes(x=Group,y=-log10(MAF))) + geom_boxplot(width=.6) + theme_pubr()

dat$Group <- factor(dat$Group,levels = c("nonOutlier","Outlier","Coloc-Outlier"))
p <- ggplot(dat,aes(x=Group,y=Beta_PCT)) + geom_violin(fill="seagreen3") + geom_boxplot(width=.1,fill="white") +
  geom_hline(yintercept = median(dat.outlier$Beta_PCT),linetype="dashed") + theme_pubr();p
