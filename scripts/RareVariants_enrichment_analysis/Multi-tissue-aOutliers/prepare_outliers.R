# combine 3UTR apa outliers and ipa outliers
# merge the tissue count matrix

setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")

dir_utr <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/output/apa/medz5"
dir_ipa <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-10-28-rareVar-IPA/output/medz_5"


library(data.table)
library(dplyr)
library(magrittr)
#load data

apa_count <- fread(paste0("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/output/apa/medz","/apa_medz.counts.txt"),header=T,sep="\t")
cat("dim(apa_count): ",dim(apa_count),"\n")
head(apa_count[1:5,1:5])

inds <- names(apa_count)[-1]

apa_outliers <- read.table(paste0(dir_utr,"/apa_medz.picked.Z_5.rm_global_chrXY.txt"),header=T,sep="\t",stringsAsFactors=F)
apa_outliers$tissue_count <- NULL

ipa_outliers <- read.table(paste0(dir_ipa,"/ipa_medz.picked.Z_5.txt"),header=T,sep="\t",stringsAsFactors=F)
cat("dim(ipa_outliers:", dim(ipa_outliers),"\n")

ipa_outliers %<>% filter(INDS %in% inds)

cat("dim(ipa_outliers:", dim(ipa_outliers),"\n")


ipa_counts <- fread(paste0(dir_ipa,"/ipa_medz.counts.txt"),header=T,sep="\t")
setdiff(inds,names(ipa_counts)[-1])
head(ipa_counts[1:5,1:5])
ipa_counts %<>% select(c("Gene",all_of(inds)))
cat("dim(ipa_counts):", dim(ipa_counts),"\n")
outlier.comb <- rbind(apa_outliers,ipa_outliers)
write.table(outlier.comb,file="./input/outliers/apa_medz.picked.Z_5.txt",quote=F,sep="\t",row.names=F)

#count.comb <- rbind(apa_count,ipa_counts)
#dim(count.comb)
#write.table(count.comb,file="./input/outliers/apa_medz.counts.txt",quote=F,sep="\t",row.names=F)

