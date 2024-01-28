# annotate CADD for all indels identfied in GTEx
library(optparse)

option_list <- list(
					make_option(c("-a","--anno"),type="character",default="",action="store",help="individual id")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(dplyr)
library(magrittr)

library(data.table)

# global settings
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/output/Count_by_Anno/1k")


# main

## load file
inds_dat <- read.table("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/gtex_inds_715.txt",
					   header=F,stringsAsFactors=F)

inds_list <- inds_dat$V1
rm(inds_dat)

for(i in 1:length(inds_list)){
	cat(inds_list[i],"\n")
	infile_1 <- paste0("./",opt$anno,"/SNP/",inds_list[i],"_counts_bygene.txt")
	infile_2 <- paste0("./",opt$anno,"/indel/",inds_list[i],"_counts_bygene.txt")

	dat1 <- fread(infile_1,header=T,sep="\t")
	names(dat1) <- c("gene_id","SNP")
	dat2 <- fread(infile_2,header=T,sep="\t")
	names(dat2) <- c("gene_id","indel")

	dat <- merge(dat1,dat2,by="gene_id")
	dat$n_variants <- dat$SNP + dat$indel
	dat %<>% select(gene_id,n_variants)

	write.table(dat,file=paste0("./",opt$anno,"/SNP_indel/",inds_list[i],"_counts_bygene.txt"),quote=F,sep="\t",row.names=F)

}

