# enrichment of RV around outliers
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue")

dir = "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/"

#library(ggplot2)
#library(cowplot)
#library(ggpubr)
options(stringsAsFactors=FALSE)

read_file = function(dir, tissue,maf,cadd,type, prefix,suffix) {
  mafstr = paste0('_MAF',maf)
  caddstr = paste0('.CADD_',cadd)
  count = read.table(paste(dir,prefix,type,mafstr,caddstr,".",tissue,suffix,sep = ""),
                     header=T,sep='\t',stringsAsFactors = F)
  count$any_variant = (count$n_variants > 0) + 0
  return(count)
}

## Function that runs proportion.ratios from enrichment.functions.R
## and also sets factor levels
counts2props = function(dir, tissues, maf, t,cadd, prefix, suffix) {
	results = data.frame(ESTIM=numeric(), CI.LOW=numeric(),CI.HIGH=numeric(),P=numeric(),TISSUE=character(),stringsAsFactors = F)

	for(tissue in tissues){
		count.df <- read_file(dir,tissue,maf,cadd,t,prefix,suffix)
		count.df <- count.df[,c('Y','any_variant')]
		summary.counts <- as.data.frame(table(count.df))
		if(nrow(summary.counts)==4 & min(summary.counts$Freq)>0){
			out.var = summary.counts$Freq[summary.counts$Y == 1 & summary.counts$any_variant == 1]
			nonout.var = summary.counts$Freq[summary.counts$Y == 0 & summary.counts$any_variant == 1]
			out.total = sum(summary.counts$Freq[summary.counts$Y == 1])
			nonout.total = sum(summary.counts$Freq[summary.counts$Y == 0])
      
			mat = c(out.var,nonout.var,out.total-out.var,nonout.total-nonout.var)
			ft <- fisher.test(matrix(mat,ncol=2))
			estimate <- ft$estimate
      # get bounds of CI
			pvalue <- ft$p.value
			max.ci <- ft$conf.int[2]
			min.ci <- ft$conf.int[1]
			dfrow = list(ESTIM=estimate, CI.LOW=min.ci, CI.HIGH=max.ci, P=pvalue, TISSUE=tissue)
			cat("ESTIM:",estimate,"CI.LOW:",min.ci,"CI.HIGH:",max.ci,"P:",pvalue,"TISSUE:",tissue,"\n")
			results <- rbind(results,dfrow)
		}else{
			dfrow = list(ESTIM=NA, CI.LOW=NA, CI.HIGH=NA, P=NA, TISSUE=tissue)
			cat("ESTIM:",estimate,"CI.LOW:",min.ci,"CI.HIGH:",max.ci,"P:",pvalue,"TISSUE:",tissue,"\n")
			results <- rbind(results,dfrow)
		}
	}
	results$TISSUE <- factor(results$TISSUE,levels = tissues)
	return(results)
}


#------------------- MAIN
## ---- different MAFs and Distance
## Set variety of variables that will be useful
MAF = '0-1'
#CADD = '0'
tissue_df <- read.table("gtex_49tissues.txt",header=F,sep="\t")
tissue_list <- as.character(tissue_df$V1)
rm(tissue_df)

TYPEs = 'SNP'
#type.names <- c("SNV","Indel")
#colors.type = c('dodgerblue3','springgreen4')
#names(colors.type) = type.names

IN.DIR = paste0(dir,'enrichment/')

# Z3,1k:CADD_0/15/25
prp_r.1k.cadd0.maf0_1 <- counts2props(IN.DIR,tissue_list,MAF,TYPEs,'0',prefix = "varCount_1k_Z5_",suffix = ".txt")
prp_r.1k.cadd15.maf0_1 <- counts2props(IN.DIR,tissue_list,MAF,TYPEs,'15',prefix = "varCount_1k_Z5_",suffix = ".txt")
prp_r.1k.cadd25.maf0_1 <- counts2props(IN.DIR,tissue_list,MAF,TYPEs,'25',prefix = "varCount_1k_Z5_",suffix = ".txt")

prp_r.1k.cadd0.maf0_1$CADD <- "CADD_0"
prp_r.1k.cadd15.maf0_1$CADD <- "CADD_15"
prp_r.1k.cadd25.maf0_1$CADD <- "CADD_25"

# Z3,2k:CADD_0/15/25

# Z3,10k:CADD_0/15/25

# combined
dat.plot <- rbind(prp_r.1k.cadd0.maf0_1,prp_r.1k.cadd15.maf0_1,prp_r.1k.cadd25.maf0_1)



saveRDS(dat.plot,file="singleTissue_aOutliers.SNV_enrichment.Z5.1k.MAF0_1.RDS")
#saveRDS(dat.plot,file="singleTissue_aOutliers.SNV_enrichment.Z4.1k.MAF_S.RDS")
#library(ggsci)

#p <- ggplot(dat.plot,aes(x=MAF,y=ESTIM,ymin=CI.LOW,ymax=CI.HIGH,color=DIST)) + 
#  geom_pointrange(linewidth=1,position = position_dodge(width = .7),size=.8) +  
#  geom_hline(yintercept = 1.0, linetype="dashed") + theme_pubr() + scale_color_npg() + 
#  facet_wrap(~Zcutoff,scales = "free_y");p
#p + scale_color_npg()
