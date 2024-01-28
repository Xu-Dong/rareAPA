# given a transcript, extract all aVariants of that transcript
library(optparse)

option_list <- list(
        make_option(c("-s","--single"),type="character",default="",action="store",help="singleton_SNPs.singletons"),
        make_option(c("-i","--inds"),type="character",default="",action="store",help="GTEX-111CU_SNPs_cadd.bed"),
        make_option(c("-o","--output"),type="character",default="",action="store",help="")
                    )

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

# load libraries
library(dplyr)
library(magrittr)
library(data.table)

# functions
paste_snp <- function(x){
	return(paste(x,collapse=":"))
}

alter_label <- function(x){
	maf <- x[5]
	label <- x[7]
	if(is.na(label)){
		return(maf)
	}else{
		return(label)
	}
}

# main

## load annotated singletons and doubletons
dat.sig <- fread(opt$single,header=T,sep="\t")
names(dat.sig) <- c("chrom","pos","label","allele","INDS")
dat.sig$pos <- as.character(dat.sig$pos)

dat.sig$snp <- apply(dat.sig[,1:2],1,paste_snp)
dat.sig %<>% select(snp,label)
cat("Dim of dat.sig:", dim(dat.sig),"\n")
## load individual variants file
dat.inds <- fread(opt$inds,header=F,sep="\t")
names(dat.inds) <- c("chrom","pos0","pos1","maf","cadd")
dat.inds$pos1 <- as.character(dat.inds$pos1)
dat.inds$snp <- apply(dat.inds[,c(1,3)],1,paste_snp)

dim(dat.inds)
cat("deduplicating...\n")
dat.inds %<>% group_by(snp) %>% mutate(min_maf=min(maf)) %>% ungroup() %>% filter(maf==min_maf)

cat("Dim of dat.inds:", dim(dat.inds),"\n")
dat.inds$min_maf <- NULL

dat.inds <- merge(dat.inds,dat.sig,by="snp",all.x=T)

cat("Dim of dat.inds after merging:", dim(dat.inds),"\n")

dat.inds$singleton <- apply(dat.inds,1,alter_label)
dat.inds %<>% distinct(.keep_all=T) %>% group_by(snp) %>% mutate(max_cadd=max(cadd)) %>% ungroup() %>% filter(cadd==max_cadd)
dat.inds$max_cadd <- NULL

dat.inds$label <- NULL
dat.inds$snp <- NULL
setorder(dat.inds,chrom,pos0)
cat("Dim of dat.inds:", dim(dat.inds),"\n")
write.table(dat.inds,file=opt$output,quote=F,sep="\t",row.names=F,col.names=F)
