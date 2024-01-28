# given a transcript plot distribution of p of leadSNP across tissues
library(optparse)

option_list <- list(
        make_option(c("-r","--ref"),type="character",default="NA",action="store",help="specify a file"),
        make_option(c("-g","--outliergenes"),type="character",default="NA",action="store",help="specify a file"),
        make_option(c("-o","--output"),type="character",default="gene_p_distr.rds",action="store",help="")
                    )

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

# load libraries
library(dplyr)
library(magrittr)
library(data.table)

extract_refid <- function(x){
	a <- strsplit(x,split="|",fixed=T)[[1]][1]
	return(strsplit(a,split=".",fixed=T)[[1]][1])
}
rm_version <- function(x){
	return(strsplit(x,split=".",fixed=T)[[1]][1])
}
# global settings
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment")


# load ref
dat.ref <- read.table(opt$ref,header=F,sep="\t")
names(dat.ref) <- c("chrom","pos0","pos1","refid","size","strand")
dat.ref$refid <- sapply(dat.ref$refid,rm_version)

# load outlier genes
dat.outlier <- read.table(opt$outliergenes,header=T,sep="\t")
dat.outlier$refid <- sapply(dat.outlier$GENE,extract_refid)

dat.ref %<>% filter(refid %in% dat.outlier$refid)

write.table(dat.ref,file=opt$output,quote=F,row.names=F,col.names=F)
