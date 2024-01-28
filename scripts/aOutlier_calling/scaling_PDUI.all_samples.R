#!/lustre/home/xdzou/src/R-4.0.5/bin/env Rscript

# R script to process peer-corrected pdui, z-transformed pduis for each gene across individuals in each tissue

# Load require packages
require(data.table)
#require(ggplot2)
#require(reshape2)


#------------- FUNCTIONS

#------------- MAIN

dir = "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project"
peer.dir = paste0(dir, '/2021-07-28-call-outliers-onlyEA/input/peer_pdui/')
#covs.dir = paste0(dir, '/data/covariates/')
#covs.dir = paste0(dir, '/data/GTEx_Analysis_v8_eQTL_covariates/') # modified by zouxd

# Read in list of tissues with eQTL data
tissue = commandArgs(trailingOnly=TRUE)

# For each tissue, read in the PDUI
# z-transform
print(tissue)
pdui <- as.data.frame(fread(paste0(peer.dir,tissue,'.pdui.peer.txt'),header=T))
rownames(pdui) <- pdui[,1]
pdui <- pdui[,-1]

pdui_ztrans = scale(t(pdui))
pdui_out <- t(pdui_ztrans)
rm(pdui_ztrans,pdui)

write.table(data.frame(Id=rownames(pdui_out),pdui_out),paste0(peer.dir, tissue, '.pdui.peer.ztrans.txt'),quote=F,sep='\t',row.names=F,col.names=T)

