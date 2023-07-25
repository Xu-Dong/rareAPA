#!/lustre/home/xdzou/src/R-4.0.5/bin/R/env Rscript

args <- commandArgs(trailingOnly=T)
if (length(args) != 2) {
  cat("Usage: R -f call_outliers_medz_peerless.R --slave --vanilla --args PDUI_tables molecular_type\n", file=stderr()) # molecular type: apa/expression
  quit(status=2)
}

dir = "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/"

# Load libraries
library(data.table)
library(reshape2)
library(plyr)
library(doMC)

###
### Setup parallel processing
###
doMC::registerDoMC(cores=8)

############ FUNCTIONS

###
### Functions for calculating number of tissues/sample and meta-analysis z-score
###
meta.n = function(values) {
	length(values) - sum(is.na(values))
}

meta.median = function(values) {
	median(values, na.rm=T)
}

meta.analysis = function(x) {
	samples=colnames(x)[3:ncol(x)]
	y = t(x[,3:ncol(x)]) # individuals (rows) x tissues (columns)
	n = apply(y, 1, meta.n)
	m1 = apply(y, 1, meta.median)
	data.frame(sample = samples, n.tissues = n, median.z = m1)
}

## Function to pick outliers from Median Z method

# For each gene we took all zscores meet the cutoff as outliers, therefore, a gene might has multiple outliers
# here the medz is already processed (remove genes with less than 5 tissues),therefore we don't need input counts
pick_outliers_v2 <- function(medz,threshold){
	outlier_idx = which(abs(medz)>threshold,arr.ind=T)
	return(data.frame(GENE=rownames(medz)[outlier_idx[,1]],
			  INDS=colnames(medz)[outlier_idx[,2]],
			  medZ=medz[outlier_idx]))
}

# For each gene we took the maximal zscore among individuals for further outlier test, this can be used for further disease gene enrichment
# here we don't specify threshold, since we will output all tested genes and their maximal zscores
pick_outliers_v3 <- function(medz,inds){
	outlier_idx = as.numeric(apply(abs(medz), 1, which.max))
	return(data.frame(INDS=factor(inds[outlier_idx], levels = inds),medZ=medz[cbind(1:nrow(medz), outlier_idx)]))
}

count_exp_gene <- function(x){
	return(length(x) - sum(is.na(x)))
}
############ MAIN

input_file <- args[1]
mol_type <- args[2] # apa/expression
###
### Loading data
###
## Load flat file with filtered and normalized expression data
data = fread(paste(dir, 'output/',input_file, sep = ''), header = T) # apa

setkey(data, Gene)

## Read in sample list
individs = colnames(data)[-c(1,2)]

## Calculate meta-analysis test statistics
# generate a data.frame with four columns(Gene,individual,tissue_count,zscore)
results = ddply(data, .(Gene), meta.analysis, .parallel = TRUE)

# filter gene ~ individual pair expressed in less than 5 tissues, make this as a parameter specified by user
library(dplyr)
tissue_threshold = 5 # change to 1
results.f <- results %>% dplyr::filter(n.tissues>=tissue_threshold)

## Unmelt results to yield data frames of tissue counts, meta Z scores, and p-values for Stouffer's method
counts = dcast(data = results, Gene ~ sample, value.var = 'n.tissues')

# write counts
write.table(counts,file=paste(dir,"output/",mol_type,"/medz_all/",mol_type,"_medz.counts.txt",sep=""),quote=F,row.names=F,sep="\t") 
rownames(counts) = counts$Gene
counts = counts[, -1]
counts = counts[, individs]

# use the filtered results.f to generate a data.frame for further outlier test
medz = dcast(data = results.f, Gene ~ sample, value.var = 'median.z')

# write medz
write.table(medz,file=paste(dir,"output/",mol_type,"/medz_all/medz_of_all_gene.",mol_type,".txt",sep=""),quote=F,row.names=F,sep="\t")
rownames(medz) <- medz$Gene
medz <- medz[,-1]

# pick outlier
Z_threshold = 3
outlier_picked = pick_outliers_v2(medz, Z_threshold)# medZ >3

# pick extreme zscore for each gene
indis = names(medz)
zscore_extreme = pick_outliers_v3(medz,indis)
zscore_extreme$GENE <- rownames(medz)
zscore_extreme <- zscore_extreme[,c(3,1,2)]

# write output
output_name = paste(dir,"output/",mol_type,"/medz_all/",mol_type,"_medz.picked.Z_",Z_threshold,".txt",sep="")
write.table(outlier_picked,file=output_name,quote=F,row.names=F,sep="\t")

#write.table(zscore_extreme,file=paste(dir,"output/",mol_type,"/medz/extreme_medz.no_filtered.",mol_type,".txt",sep=""),quote=F,row.names=F,sep="\t")

