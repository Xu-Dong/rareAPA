#!/lustre/home/xdzou/src/R-4.0.5/bin/env Rscript

# extract covariates for only European individuals for each tissue

# usage: Rscript extract_european_inds.covariates.R tissue.v8.covariates.txt


basedir <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA"
# load args
args <- commandArgs(trailingOnly=TRUE)

dat <- read.table(args[1],header=T,sep="\t",stringsAsFactors=F)
inds <- names(dat)[-1]
inds <- gsub("\\.","-",inds)
names(dat) <- c("ID",inds)
# load EA list
eaList <- read.table("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/input/individuals.onlyEA.txt", header=F, stringsAsFactors=F)

names(eaList) <- c("indis")

inds.ov <- intersect(eaList$indis,inds)

# output
dat.out <- dat[,c("ID",inds.ov)]

write.table(dat.out,file=paste0(basedir,"/input/GTExV8_COV_EA/",basename(args[1])),quote=F,sep="\t",row.names=F)


