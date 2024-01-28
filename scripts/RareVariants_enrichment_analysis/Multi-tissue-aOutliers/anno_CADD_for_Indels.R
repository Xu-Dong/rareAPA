# annotate CADD for all indels identfied in GTEx
library(optparse)

option_list <- list(
					make_option(c("-c","--chrom"),type="integer",default="chr1",action="store",help="chromosome")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(dplyr)
library(magrittr)

library(data.table)

setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")
# functions
extract_allele <- function(x){
	return(strsplit(x,split=":",fixed=T)[[1]][1])
}

paste_indel <- function(x){
	return(paste(x,collapse=":"))
}

paste_chr_to_pos <- function(x){
	return(paste0("chr",x))
}

# load all gtex indels

chromosome <- paste0("chr",opt$chrom)

gtex.indel <- fread("./input/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze_indels.frq",header=T,sep="\t")
gtex.indel <- gtex.indel[,c(1,2,5,6)]
names(gtex.indel) <- c("chrom","pos","ref","alt")

gtex.indel$ref <- sapply(gtex.indel$ref,extract_allele)
gtex.indel$alt <- sapply(gtex.indel$alt,extract_allele)

gtex.indel$pos <- as.character(gtex.indel$pos)

gtex.indel %<>% filter(chrom==chromosome)

all_gtex_indels <- apply(gtex.indel,1,paste_indel)

rm(gtex.indel)


# load indel cadd (by chromosome)
cadd.indel <- fread(paste0("./input/tmp/InDel_cadd.",opt$chrom,".tsv"),header=F,sep="\t")
names(cadd.indel) <- c("chrom","pos","ref","alt","CADD")

cadd.indel$chrom <- sapply(cadd.indel$chrom,paste_chr_to_pos)
cadd.indel$pos <- as.character(cadd.indel$pos)

cadd.indel$indel <- apply(cadd.indel[,1:4],1,paste_indel)

cadd.indel %<>% filter(indel %in% all_gtex_indels) %>% select(chrom,pos,CADD)

write.table(cadd.indel,file=paste0("./input/tmp/GTEx_indel_cadd.",opt$chrom,".txt"),quote=F,sep="\t",row.names=F)
