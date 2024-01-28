# combine 3UTR apa outliers and ipa outliers
# merge the tissue count matrix

setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")

dir <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/output/GTEx_v8_outlier_calls"


library(data.table)
library(dplyr)
library(magrittr)


rm_version <- function(x){
	return(strsplit(x,split=".",fixed=T)[[1]][1])
}

#load data

eoutlier <- read.table(paste0(dir,"/eoutliers_medz.picked.Z_3.txt"),header=T,sep="\t")
soutlier <- read.table(paste0(dir,"/soutliers_medz.picked.Z_3.txt"),header=T,sep="\t")

eoutlier$GENE <- sapply(eoutlier$GENE,rm_version)
soutlier$GENE <- sapply(soutlier$GENE,rm_version)

#write.table(eoutlier,file="./input/eOutliers_medz.picked.Z_3.txt",quote=F,sep="\t",row.names=F)
#write.table(soutlier,file="./input/sOutliers_medz.picked.Z_3.txt",quote=F,sep="\t",row.names=F)


union_gene_set <- unique(c(eoutlier$GENE,soutlier$GENE))
inds_dat <- read.table(paste0(dir,"/individuals.onlyEA.txt"),header=F,stringsAsFactors=F)
inds_list <- inds_dat$V1
rm(inds_dat)
N_gene <- length(union_gene_set)
N_inds <- length(inds_list)
dat.m <- matrix(rep(10,times=N_gene*N_inds),nrow=N_gene,ncol=N_inds)

dat.m <- as.data.frame(cbind(union_gene_set,dat.m))
colnames(dat.m) <- c("Gene",inds_list)

#write.table(dat.m,file="./input/exp_medz.counts.txt",quote=F,sep="\t",row.names=F)


# load genbody region
cat("Number of union genes:", N_gene,"\n")

dat.gb <- fread("./input/gencode.v26.genes.v8.patched_contigs.bed",header=F,sep="\t")
dat.gb$V4 <- sapply(dat.gb$V4,rm_version)

x <- intersect(dat.gb$V4,union_gene_set)
cat("Number of overlapped genes:", length(x),"\n")

dat.gb %<>% filter(V4 %in% union_gene_set)
dim(dat.gb)
write.table(dat.gb,file="./input/eOutlier_genes_gb.txt",quote=F,sep="\t",row.names=F,col.names=F)


