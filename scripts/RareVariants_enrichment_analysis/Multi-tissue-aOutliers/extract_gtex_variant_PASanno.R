# annotate CADD for all indels identfied in GTEx
library(optparse)

option_list <- list(
					make_option(c("-i","--inds"),type="character",default="",action="store",help="individual id"),
					make_option(c("-v","--vartype"),type="character",default="SNP",action="store",help="SNP/indel")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(dplyr)
library(magrittr)

library(data.table)

# global settings
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")
dir <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-30-RV-Annotation-VEP/output/individual_anno_bySite/"


# main

## load file
infile <- paste0(dir,opt$vartype,"/",opt$inds,".anno_",opt$vartype,".bed")
dat <- fread(infile,header=F,sep="\t")

## trim
dat <- dat[,c(1,2,3,4,19)]
names(dat) <- c("chrom","pos0","pos1","maf","PAS")

dim(dat)
## removing common variants
dat %<>% filter(maf<0.01)
dim(dat)


write.table(dat,file=paste0("./input/GTEx_Variants_with_PASanno/",opt$inds,"_",opt$vartype,"_anno.bed"),quote=F,sep="\t",row.names=F,col.names=F)
