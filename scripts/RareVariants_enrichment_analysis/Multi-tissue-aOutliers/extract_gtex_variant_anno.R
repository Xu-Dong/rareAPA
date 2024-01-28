# annotate CADD for all indels identfied in GTEx
library(optparse)

option_list <- list(
					make_option(c("-i","--inds"),type="character",default="",action="store",help="individual id"),
					make_option(c("-v","--vartype"),type="character",default="SNP",action="store",help="SNP/indel")
					)

opt <- parse_args(OptionParser(option_list=option_list,usage="usage: %prog [options]"))

library(dplyr)
library(magrittr)

library(data.table)

# global settings
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")
dir <- "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-30-RV-Annotation-VEP/output/individual_anno_bySite/"


# main

## load file
infile <- paste0(dir,opt$vartype,"/",opt$inds,".anno_",opt$vartype,".bed")
dat <- fread(infile,header=F,sep="\t")

## trim
dat <- dat[,c(1,2,3,4,7,8,9,13,14,15,16,17)]
names(dat) <- c("chrom","pos0","pos1","maf","nTFBS","nPromoter","nEnhancer","mamPhastCons","mamPhyloP","GerpN","GerpS","VEP")

dim(dat)
## removing common variants
dat %<>% filter(maf<0.01)
dim(dat)

dat$VEP_SO <- NA
dat$VEP_SO[dat$VEP=="transcript_ablation"] <- "Transcript_ablation/amplification"
dat$VEP_SO[dat$VEP=="splice_acceptor_variant"] <- "Splice"
dat$VEP_SO[dat$VEP=="splice_donor_variant"] <- "Splice"
dat$VEP_SO[dat$VEP=="splice_region_variant"] <- "Splice"
dat$VEP_SO[dat$VEP=="stop_gained"] <- "Stop"
dat$VEP_SO[dat$VEP=="stop_lost"] <- "Stop"
dat$VEP_SO[dat$VEP=="frameshift_variant"] <- "Frame_shift"
dat$VEP_SO[dat$VEP=="start_lost"] <- "Stop"
dat$VEP_SO[dat$VEP=="transcript_amplification"] <- "Transcript_ablation/amplification"
dat$VEP_SO[dat$VEP=="inframe_insertion"] <- "Coding"
dat$VEP_SO[dat$VEP=="inframe_deletion"] <- "Coding"
dat$VEP_SO[dat$VEP=="missense_variant"] <- "Coding"
dat$VEP_SO[dat$VEP=="incomplete_terminal_codon_variant"] <- "Coding"
dat$VEP_SO[dat$VEP=="start_retained_variant"] <-"Coding" 
dat$VEP_SO[dat$VEP=="stop_retained_variant"] <- "Coding"
dat$VEP_SO[dat$VEP=="synonymous_variant"] <- "otherCoding"
dat$VEP_SO[dat$VEP=="coding_sequence_variant"] <- "otherCoding"
dat$VEP_SO[dat$VEP=="mature_miRNA_variant"] <- "Regulatory-Noncoding" # new modified,need re-run
dat$VEP_SO[dat$VEP=="5_prime_UTR_variant"] <- "UTR_5"
dat$VEP_SO[dat$VEP=="3_prime_UTR_variant"] <- "UTR_3"
dat$VEP_SO[dat$VEP=="non_coding_transcript_exon_variant"] <- "otherNoncoding"
dat$VEP_SO[dat$VEP=="intron_variant"] <- "Intron"
dat$VEP_SO[dat$VEP=="NMD_transcript_variant"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="non_coding_transcript_variant"] <- "otherNoncoding"
dat$VEP_SO[dat$VEP=="upstream_gene_variant"] <- "otherNoncoding"
dat$VEP_SO[dat$VEP=="downstream_gene_variant"] <- "otherNoncoding"
dat$VEP_SO[dat$VEP=="TFBS_ablation"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="TFBS_amplification"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="TF_binding_site_variant"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="regulatory_region_ablation"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="regulatory_region_amplification"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="feature_elongation"] <- "Others"
dat$VEP_SO[dat$VEP=="regulatory_region_variant"] <- "Regulatory-Noncoding"
dat$VEP_SO[dat$VEP=="feature_truncation"] <- "Others"
dat$VEP_SO[dat$VEP=="intergenic_variant"] <- "Intergenic"


write.table(dat,file=paste0("./input/GTEx_Variants_with_Anno/",opt$inds,"_",opt$vartype,"_anno.bed"),quote=F,sep="\t",row.names=F,col.names=F)
