# annotate CADD for all indels identfied in GTEx
library(optparse)

option_list <- list(
					make_option(c("-i","--input"),type="character",default="",action="store",help="individual id"),
					make_option(c("-t","--tissue"),type="character",default="",action="store",help="SNP/indel")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(dplyr)
library(magrittr)

library(data.table)

# global settings
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")
dir <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/"

# main

## load file
cat("loading the big matrix...\n")
dat <- fread(opt$input,header=T,sep="\t")
cat("Finally, loading completed!\n")
tissues <- read.table(opt$tissue,header=F)
tissue_list <- tissues$V1
for(i in 1:length(tissue_list)){
	cat(tissue_list[i],'\n')
	dat.tissue <- dat %>% filter(Tissue==tissue_list[i])
	dat.tissue$Tissue <- NULL
	write.table(dat.tissue,file=paste0(dir,"single_tissue/Zmatrix/",tissue_list[i],".z_mat.txt"),quote=F,sep="\t",row.names=F)
}

