#!/bin/Rscript

library(dplyr)
library(magrittr)

library(CMplot)
library(ggplot2)
library(ggpubr)

manhattan.df <- readRDS("manhattan_plot.aOutliers.RDS")
outliers <- manhattan.df %>% filter(abs(Z)>3)
dat.p1 <- manhattan.df%>% select(SNP,chrom,pos,Z)
names(dat.p1) <- c("ID","CHR","Pos_TSS","Z")
color_code <- data.frame(Class=c("3UTR","Intronic"),Color_code=c("#3B4992FF","#EE0000FF"))
outliers <- merge(outliers,color_code,by="Class",all.x = T)

CMplot(dat.p1,plot.type = "m",
       multracks=FALSE,
       band = 1,
       LOG10 = FALSE,
       ylim = c(-15,20),
       ylab="Median Z-score",
       amplify = FALSE,
       width=14,
       height = 6,
       col=c("grey30","grey60"),
       threshold=c(-3,3),
       threshold.lty = c(2,2),
       threshold.lwd = c(1,1),
       highlight = outliers$SNP,
       highlight.col = outliers$Color_code,
       file="pdf",
       file.output=TRUE,
       cex = 1.0,
       pch = 20,
       verbose=TRUE
)
