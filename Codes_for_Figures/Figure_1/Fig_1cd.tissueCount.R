# distribution of tissue count of each aOutlier (3'UTR aOutlier vs. Intronic aOutlier
rm(list=ls())

library(dplyr)
library(data.table)
library(magrittr)

# load multi-tissue aOutliers
dat <- fread("aOutliers_multi-tissue.txt",header=T,sep="\t") # multi-tissue aOutliers is available through our website:http://bioinfo.szbl.ac.cn/rareAPA/Download.php
names(dat) <- c("Class","Transcript","Location","Symbol","INDS","tissue_count","medZ")
dat.3utr <- dat %>% filter(Class=="3'UTR")
dat.intron <- dat %>% filter(Class=="Intronic")

# - barplot of tissue count
df.tissue.apa <- as.data.frame(table(dat.3utr$tissue_count))
df.tissue.ipa <- as.data.frame(table(dat.intron$tissue_count))
names(df.tissue.apa) <- c("TissueCount","Freq_apa")
names(df.tissue.ipa) <- c("TissueCount","Freq_ipa")
df.tissue.apa$Prop_apa <- df.tissue.apa$Freq_apa/sum(df.tissue.apa$Freq_apa)
df.tissue.ipa$Prop_ipa <- df.tissue.ipa$Freq_ipa/sum(df.tissue.ipa$Freq_ipa)
df.plot <- data.frame(TissueCount=5:36)
df.plot <- merge(df.plot,df.tissue.apa[,c(1,3)],by="TissueCount",all.x = T)
df.plot <- merge(df.plot,df.tissue.ipa[,c(1,3)],by="TissueCount",all.x = T)
df.plot$Prop_apa[is.na(df.plot$Prop_apa)] <- 0
df.plot$Prop_ipa[is.na(df.plot$Prop_ipa)] <- 0
df.plot$TissueCount <- factor(df.plot$TissueCount,levels = 5:36)

library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
df.plot <- melt(df.plot,id.vars = "TissueCount")
p <- ggplot(df.plot,aes(x=TissueCount,y=value,fill=variable)) + geom_bar(stat = "identity",width = .5,position = "dodge") +
  theme_pubr() + scale_fill_aaas();p

# - boxplot of MeanTissue
tissue_stat.apa <- dat.apa %>% group_by(GENE) %>% summarise(meanTissue=mean(tissue_count))
tissue_stat.ipa <- dat.ipa %>% group_by(GENE) %>% summarise(meanTissue=mean(tissue_count))
plot.df <- data.frame(Count=c(tissue_stat.apa$meanTissue,tissue_stat.ipa$meanTissue),
                      Group=c(rep("APA",dim(tissue_stat.apa)[1]),rep("IPA",dim(tissue_stat.ipa)[1])))
p.vl <- ggplot(plot.df,aes(x=Group,y=Count)) + geom_boxplot(aes(fill=Group),width=.5) + 
  theme_pubr() + scale_y_continuous(breaks=c(5,10,15,20,25,30,35,40),labels=c(5,10,15,20,25,30,35,40)) +
  xlab("") + ylab("Acrossing tissue number")

p.vl.aaas <- p.vl + scale_fill_aaas();p.vl.aaas

# test the difference
wilcox.test(Count~Group,data = plot.df,alternative="less")
