#check global medz outliers

setwd("C:/WorkSpace/ScienceData/Projects/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/")


# -- load libraries
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
# -- functions
count_nonNA <- function(x){
  return(sum(!is.na(x)))
}

figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=14, hjust=0.5), 
               text = element_text(size=14),axis.text=element_text(size=14),
               axis.title = element_text(size=16)))
}
# -- load data
# apa
aoutliers <- read.table("./output/apa/medz/apa_medz.picked.Z_3.txt",header=TRUE, sep = "\t", 
                        stringsAsFactors = FALSE)
medz.apa <- fread("./output/apa/medz/medz_of_all_gene.apa.txt",header=TRUE, sep="\t")
# expression
eoutliers <- read.table("./output/expression/medz/expression_medz.picked.Z_3.txt",stringsAsFactors = F,
                        header = T, sep = "\t")
medz.exp <- fread("./output/expression/medz/medz_of_all_gene.expression.txt",header=TRUE, sep = "\t")

# -- count tested genes
# apa
count_tested_genes.apa <- apply(medz.apa[,2:dim(medz.apa)[2]],2,count_nonNA)
names(count_tested_genes.apa) <- colnames(medz.apa)[-1]
count_tested_genes.apa.df <- data.frame(INDS=names(count_tested_genes.apa),
                                   count=count_tested_genes.apa)

# expression
count_tested_genes.exp <- apply(medz.exp[,2:dim(medz.exp)[2]],2,count_nonNA)
names(count_tested_genes.exp) <- colnames(medz.exp)[-1]
count_tested_genes.exp.df <- data.frame(INDS=names(count_tested_genes.exp),
                                        count=count_tested_genes.exp)

# -- count outlier genes per individual and filter global outliers
# apa
count.df.apa <- data.frame(table(aoutliers$INDS))
names(count.df.apa) <- c("INDS","Freq")

count.df.apa <- merge(count.df.apa,count_tested_genes.apa.df,by="INDS")
count.df.apa$Prop <- -log10(count.df.apa$Freq/count.df.apa$count)

# plot box
boxplot.apa <- ggplot(count.df.apa,aes(x="",y=Prop)) + geom_boxplot() + theme_pubr()
boxplot.apa <- boxplot.apa + xlab("") + ylab("log10 (Proportion)\n") + 
  ggtitle(label = "Proportion of apaOutlier genes per individual") + 
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=15),
        plot.title = element_text(size=16,face = "bold",hjust = .5))

q1 <- quantile(count.df.apa$Prop,0.25)
q3 <- quantile(count.df.apa$Prop,0.75)
iqr <- q3 - q1

exclude.inds.apa <- count.df.apa %>% filter(Prop>q3+1.5*iqr | Prop<q1-1.5*iqr) # 11 individuals are global outliers

aoutliers.f <- aoutliers %>% filter(!INDS %in% exclude.inds.apa$INDS)

write.table(aoutliers.f,file = "./output/apa/medz/apa_medz.picked.Z_3.rm_global.txt",
            quote = F, row.names = F, sep = "\t")

count.df.apa$global_inds <- "NO"
count.df.apa$global_inds[count.df.apa$INDS %in% exclude.inds.apa$INDS] <- "YES"
count.df.apa %<>% select(INDS,Freq,global_inds)
write.table(count.df.apa,file = "./output/apa/medz/apa_medz.picked_count_per_inds.txt",
            quote=F,row.names = F, sep = "\t")

# gene expression 
count.df.exp <- data.frame(table(eoutliers$INDS))
names(count.df.exp) <- c("INDS","Freq")

count.df.exp <- merge(count.df.exp,count_tested_genes.exp.df,by="INDS")
count.df.exp$Prop <- -log10(count.df.exp$Freq/count.df.exp$count)

# plot box
boxplot.exp <- ggplot(count.df.exp,aes(x="",y=Prop)) + geom_boxplot() + theme_pubr()
boxplot.exp <- boxplot.exp + xlab("") + ylab("log10 (Proportion)\n") + 
  ggtitle(label = "Proportion of eOutlier genes per individual") + 
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=15),
        plot.title = element_text(size=16,face = "bold",hjust = .5))

q1 <- quantile(count.df.exp$Prop,0.25)
q3 <- quantile(count.df.exp$Prop,0.75)
iqr <- q3 - q1

exclude.inds.exp <- count.df.exp %>% filter(Prop>q3+1.5*iqr | Prop<q1-1.5*iqr) # 19 individuals are global outliers

eoutliers.f <- eoutliers %>% filter(!INDS %in% exclude.inds.exp$INDS)

write.table(aoutliers.f,file = "./output/expression/medz/exp_medz.picked.Z_3.rm_global.txt",
            quote = F, row.names = F, sep = "\t")

count.df.exp$global_inds <- "NO"
count.df.exp$global_inds[count.df.exp$INDS %in% exclude.inds.exp$INDS] <- "YES"
count.df.exp %<>% select(INDS,Freq,global_inds)
write.table(count.df.exp,file = "./output/expression/medz/exp_medz.picked_count_per_inds.txt",
            quote=F,row.names = F, sep = "\t")


# -- combine plots
library(cowplot)

cowplot::plot_grid(boxplot.apa,boxplot.exp,ncol = 2,labels = c("A","B"))
