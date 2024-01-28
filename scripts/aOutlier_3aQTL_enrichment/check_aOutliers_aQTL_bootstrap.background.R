library(optparse)

option_list <- list(
					make_option(c("-t","--tissue"),type="character",default="NA",action="store",help="specify a gene")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))


setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-01-15-aQTL_aOutlier")

library(data.table)
library(dplyr)
library(magrittr)


# load multi-tissue
#df_mult <- fread("./input/apa_medz.picked.Z_2.rm_global_chrX_Y.txt",header=T,sep="\t")
df_mult <- fread("./input/apa_medz.picked.Z_3.rm_global_and_chrX.add_tissueCount.txt",header=T,sep="\t")
#df_mult <- fread("./input/apa_medz.picked.Z_4.rm_global_chrXY.txt",header=T,sep="\t")
#df_mult <- fread("./input/apa_medz.picked.Z_5.rm_global_chrXY.txt",header=T,sep="\t")

df_all <- fread("./input/All_tested_genes.txt",header=F,sep="\t")

m_aoutlierG <- unique(df_mult$GENE)
all_uniqG <- unique(df_all$V1)

# load single z 3

df_s <- fread("./input/aOutliers_singlez_picked.Z3.txt",header=T,sep="\t")

df_s %<>% filter(GENE %in% m_aoutlierG)


# load aQTL profile

df_aqtl <- fread("./input/aGene_profile.FDR_0.05.txt",header=T,sep="\t")

curr_tissue <- opt$tissue
aGene_list <- df_aqtl$gene[!is.na(df_aqtl[[curr_tissue]])]
if(curr_tissue=="Brain_Spinal_cord_cervical_c-1"){
	curr_tissue <- "Brain_Spinal_cord_cervical_c1"
}
aOutlier_genes <- unique(df_s$GENE[df_s$TISSUE==curr_tissue])
N <- length(aOutlier_genes)
cat(curr_tissue,"\t", N,"\n")
set.seed(1000)
ratio_list <- c()
for(i in 1:1000){
	aOutlier_sub <- sample(all_uniqG,N,replace=TRUE)
	m <- length(intersect(aOutlier_sub,aGene_list))
	ratio_list[i] <- m/N
}
Ratio_Mean <- mean(ratio_list,na.rm=T)
Lower.CI <- quantile(ratio_list,0.05)
Upper.CI <- quantile(ratio_list,0.95)

df_out <- data.frame(Tissue=c(curr_tissue), Mean=c(Ratio_Mean),Lower_CI = c(Lower.CI),Upper_CI = c(Upper.CI))

write.table(df_out,file=paste0("./output/bootstrap_z3_bg/",curr_tissue,".bootstrap_ratio.Z3.txt"),quote=F,sep="\t",row.names=F)
