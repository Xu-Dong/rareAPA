# compare GC/AT content in 3'UTR between aOutlier genes and eOutlier genes
library(data.table)
library(dplyr)
library(magrittr)


gc_dat <- fread("data/02.GC_content.txt",header=T,sep = "\t")
aoutlier_dat <- fread("data/apa_medz.picked.Z_3.rm_global_and_chrX.add_tissueCount.txt",header=T,sep = "\t")
eoutlier_dat <- fread("src/eoutliers_medz.picked.Z_3.txt",header=T,sep = "\t")
id_con <- fread("src/gencode.v26.geneid_to_genename.txt",header=F,sep = "\t")

aoutlier_dat$ID <-lapply(aoutlier_dat$GENE, function(x){strsplit(x,split="|",fixed=T)[[1]][1]})
aoutlier_dat$ID <- as.character(aoutlier_dat$ID)
aoutlier_dat_gc <- merge(aoutlier_dat,gc_dat,by="ID",all.x = TRUE)
aoutlier_dat_gc <- aoutlier_dat_gc[,-c(3:5)]


eoutlier_dat_con <- merge(eoutlier_dat,id_con,by.x="GENE",by.y="V1",all.x = TRUE)
eoutlier_dat_gc <- merge(eoutlier_dat_con,gc_dat,by.x="V2",by.y="Gene",all.x = TRUE)
eoutlier_dat_gc <- eoutlier_dat_gc[,-c(3:5)]

aoutlier_dat_gc %<>% distinct(.keep_all = T)
aoutlier_dat_gcscore <- aoutlier_dat_gc
aoutlier_dat_gcscore %<>% group_by(Gene) %>% mutate(max_GCscore=max(GCscore)) %>% ungroup() %>% filter(GCscore==max_GCscore)

eoutlier_dat_gc %<>% distinct(.keep_all = T)
eoutlier_dat_gcscore <- eoutlier_dat_gc
eoutlier_dat_gcscore %<>% group_by(V2) %>% mutate(max_GCscore=max(GCscore)) %>% ungroup() %>% filter(GCscore==max_GCscore)

### analysis from GC ratio
aoutlier_dat_gc$ratio <- aoutlier_dat_gc$GCscore/aoutlier_dat_gc$UTR_length
aoutlier_dat_ratio <- aoutlier_dat_gc
aoutlier_dat_ratio %<>% group_by(Gene) %>% mutate(max_ratio=max(ratio)) %>% ungroup() %>% filter(ratio==max_ratio)

eoutlier_dat_gc$ratio <- eoutlier_dat_gc$GCscore/eoutlier_dat_gc$UTR_length
eoutlier_dat_ratio <- eoutlier_dat_gc
eoutlier_dat_ratio %<>% group_by(V2) %>% mutate(max_ratio=max(ratio)) %>% ungroup() %>% filter(ratio==max_ratio)

aoutlier_dat_ratio$at_ratio <- aoutlier_dat_ratio$ATscore/aoutlier_dat_ratio$UTR_length
eoutlier_dat_ratio$at_ratio <- eoutlier_dat_ratio$ATscore/eoutlier_dat_ratio$UTR_length

aoutlier_dat_ratio$Group <- "aoutlier"
eoutlier_dat_ratio$Group <- "eoutlier"
names(eoutlier_dat_ratio)[1] <- "Gene"

dat.m <- rbind(aoutlier_dat_ratio,eoutlier_dat_ratio)
dat.m <- as.data.frame(dat.m)
dat.m$Group <- factor(dat.m$Group,levels = c("eoutlier","aoutlier"))

library(cowplot)
library (ggplot2)
library(ggsci)
library(ggpubr)
#pdf("aOutlier_vs_eOutlier_GCscore_ratio.color.pdf",width=3.4, height=5.6)
p1 <- ggplot(dat.m,aes(x=Group,y=ratio,fill=Group)) + geom_boxplot(width=.5) + theme_pubr() + 
  scale_fill_manual(breaks = c("eoutlier","aoutlier"),values = c("#4DBBD5FF","#E64B35FF"))

p2<- ggplot(dat.m,aes(x=Group,y=at_ratio,fill=Group)) + geom_boxplot(width=.5) + theme_pubr() + 
  scale_fill_manual(breaks = c("eoutlier","aoutlier"),values = c("#4DBBD5FF","#E64B35FF")) +
  ylab("ratio")
#dev.off()

#save.image(file="project_image.RData")
wilcox.test(aoutlier_dat_ratio$ratio,eoutlier_dat_ratio$ratio,alternative = "less")
wilcox.test(aoutlier_dat_ratio$at_ratio,eoutlier_dat_ratio$at_ratio,alternative = "greater")
