library(dplyr)
library(magrittr)
library(ggplot2)

library(ggpubr)
library(cowplot)
library(RColorBrewer)
# load raw medZ matrix 
medZ <- readRDS("./example_data/medZ_matrix.RDS")

# 11 individuals are global outleirs (need to be removed later)
gb_inds <- c("GTEX-11EQ8","GTEX-132AR","GTEX-13IVO","GTEX-13NZ8","GTEX-14E6C","GTEX-1S5ZU","GTEX-1S82Z","GTEX-P44H","GTEX-POYW","GTEX-QMR6","GTEX-QXCU")

# ------- plot for single aOutleir gene ----------------
gene <- "NM_172231.4|SUGP1|chr19|-"

gene.z <- as.numeric(medZ[medZ$Gene==gene,-1])
names(gene.z) <- names(medZ)[-1]

gene.z <- gene.z[!names(gene.z) %in% gb_inds]
genez.df <- data.frame(INDS=names(gene.z),Zscore=gene.z)
genez.df %<>% filter(!is.na(Zscore))

genez.df$isOutlier <- "NO"
genez.df$isOutlier[abs(genez.df$Zscore)>3] <- "YES"
genez.df$rank_z <- rank(genez.df$Zscore,ties.method = "random")
genez.df$sign <- ifelse(genez.df$isOutlier=="YES",genez.df$INDS,NA)
library(ggrepel)

genename <- strsplit(gene,split = "|",fixed = T)[[1]][2]
genez.df <- genez.df %>% filter(!is.na(Zscore))
p <- ggplot(genez.df,aes(x=rank_z,y=Zscore,color=isOutlier)) + 
  geom_point(alpha=.8, size=2.5) + theme_pubr() + 
  xlab("Samples ranked by normalized PDUI") + ylab("Normalized PDUI") + 
  scale_color_manual(name="",values = c("gray","red")) + 
  theme(legend.position = "") + 
  geom_label_repel(aes(label=sign), color="black", 
                   box.padding=unit(0.5, "lines"), point.padding=unit(0.8, "lines"), 
                   segment.colour = "grey50")
