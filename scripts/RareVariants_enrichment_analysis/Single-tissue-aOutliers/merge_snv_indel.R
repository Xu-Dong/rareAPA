# annotate CADD for all indels identfied in GTEx
library(optparse)

option_list <- list(
					make_option(c("-t","--tissue"),type="character",default="",action="store",help="SNP/indel")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(dplyr)
library(magrittr)

library(data.table)

# global settings
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue")

# main

## load file
cat("loading the big matrix...\n")
df_snv <- fread(paste0("./enrichment/varCount_1k_Z3_SNPs_MAF0-1.CADD_15.",opt$tissue,".txt"),header=T,sep="\t")
df_indel <- fread(paste0("./enrichment/varCount_1k_Z3_indels_MAF0-1.CADD_15.",opt$tissue,".txt"),header=T,sep="\t")
cat("Finally, loading completed!\n")

x <- df_snv$n_variants
y <- df_indel$n_variants

df_snv$n_variants <- x+y

rm(x,y)

write.table(df_snv,file=paste0("./enrichment_snv_indel/varCount_1k_Z3_RVs_MAF0-1.CADD_15.",opt$tissue,".txt"),quote=F,sep="\t",row.names=F)
