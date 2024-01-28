#check global medz outliers
library(optparse)
option_list <- list(
					make_option(c("-d","--dir"),type="character", action="store",help="filename"),
					make_option(c("-z","--zscore"),type="character",action="store",help="filename")) 

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

# change to the base dir
setwd(opt$dir)


# -- load libraries
library(dplyr)
library(data.table)
library(magrittr)
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
aoutliers <- read.table(paste0("./output/apa/medz_all/apa_medz.picked.Z_3.txt"),header=TRUE, sep = "\t",stringsAsFactors = FALSE)
medz.apa <- fread(paste0("./output/apa/medz_all/medz_of_all_gene.apa.txt"),header=TRUE, sep="\t")

# -- count tested genes
# apa
count_tested_genes.apa <- apply(medz.apa[,2:dim(medz.apa)[2]],2,count_nonNA)
names(count_tested_genes.apa) <- colnames(medz.apa)[-1]
count_tested_genes.apa.df <- data.frame(INDS=names(count_tested_genes.apa),
                                   count=count_tested_genes.apa)


# -- count outlier genes per individual and filter global outliers
# apa
count.df.apa <- data.frame(table(aoutliers$INDS))
names(count.df.apa) <- c("INDS","Freq")

count.df.apa <- merge(count.df.apa,count_tested_genes.apa.df,by="INDS")
count.df.apa$Prop <- -log10(count.df.apa$Freq/count.df.apa$count)
#saveRDS(count.df.apa,file=paste0("./output/apa/medz1_5/count.df.apa.RDS"))
# plot box
#library(ggplot2)
#library(ggpubr)
#pdf(paste0("./output/recall_ipaOutliers/medz_",opt$z,"/boxplot_global_outliers.Z_",opt$z,".pdf"),width=4,height=4)
#boxplot.apa <- ggplot(count.df.apa,aes(x="",y=Prop)) + geom_boxplot() + theme_pubr()
#boxplot.apa <- boxplot.apa + xlab("") + ylab("log10 (Proportion)\n") + 
#  ggtitle(label = "Proportion of ipaOutlier genes per individual") + 
#  theme(axis.title = element_text(size=16),
#        axis.text = element_text(size=15),
#        plot.title = element_text(size=16,face = "bold",hjust = .5))
#
#print(boxplot.apa)
#dev.off()

q1 <- quantile(count.df.apa$Prop,0.25)
q3 <- quantile(count.df.apa$Prop,0.75)
iqr <- q3 - q1

exclude.inds.apa <- count.df.apa %>% filter(Prop>q3+1.5*iqr | Prop<q1-1.5*iqr) # 11 individuals are global outliers
write.table(exclude.inds.apa,file="./output/apa/medz_all/global_outliers.txt",quote=F,sep="\t",row.names=F)
aoutliers.f <- aoutliers %>% filter(!INDS %in% exclude.inds.apa$INDS)

write.table(aoutliers.f,file = "./output/apa/medz_all/apa_medz.picked.Z_3.rm_global.txt", quote = F, row.names = F, sep = "\t")

count.df.apa$global_inds <- "NO"
count.df.apa$global_inds[count.df.apa$INDS %in% exclude.inds.apa$INDS] <- "YES"
count.df.apa %<>% select(INDS,Freq,global_inds)
write.table(count.df.apa,file = paste0("./output/apa/medz_all/apa_medz.picked_count_per_inds.txt"),quote=F,row.names = F, sep = "\t")
