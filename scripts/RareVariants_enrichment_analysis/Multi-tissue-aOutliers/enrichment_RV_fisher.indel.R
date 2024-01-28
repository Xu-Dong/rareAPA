# enrichment of RV around outliers
setwd("/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier")

dir = "/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/"

#library(ggplot2)
#library(cowplot)
#library(ggpubr)
options(stringsAsFactors=FALSE)

read_file = function(dir, maf,type, prefix,suffix) {
  mafstr = paste0('_MAF',maf)
  count = read.table(paste(dir,prefix,type,mafstr,suffix,sep = ""),
                     header=T,sep='\t',stringsAsFactors = F)
  count$any_variant = (count$n_variants > 0) + 0
  return(count)
}

## Function that runs proportion.ratios from enrichment.functions.R
## and also sets factor levels
counts2props = function(dir, mafs, t, prefix, suffix) {
	results = data.frame(ESTIM=numeric(), CI.LOW=numeric(),CI.HIGH=numeric(),P=numeric(),MAF=character(),stringsAsFactors = F)

	for(m in mafs){
		count.df <- read_file(dir,m,t,prefix,suffix)
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
			dfrow = list(ESTIM=estimate, CI.LOW=min.ci, CI.HIGH=max.ci, P=pvalue, MAF=m)
			cat("ESTIM:",estimate,"CI.LOW:",min.ci,"CI.HIGH:",max.ci,"P:",pvalue,"MAF:",m,"\n")
			results <- rbind(results,dfrow)
		}else{
			dfrow = list(ESTIM=NA, CI.LOW=NA, CI.HIGH=NA, P=NA, MAF=m)
			cat("ESTIM:",NA,"CI.LOW:",NA,"CI.HIGH:",NA,"P:",NA,"MAF:",m,"\n")
			results <- rbind(results,dfrow)
		}
	}
	results$MAF <- factor(results$MAF,levels = mafs)
	return(results)
}


#------------------- MAIN
## ---- different MAFs and Distance
## Set variety of variables that will be useful
MAFs = c('-S','0-1', '1-5', '5-10', '10-25')

TYPEs = 'indels'
#type.names <- c("SNV","Indel")
#colors.type = c('dodgerblue3','springgreen4')
#names(colors.type) = type.names

IN.DIR = paste0(dir,'output/enrichment/')

# Z3,1k:CADD_0/15/25
prp_r.1k.cadd0 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_1k_Z3_",suffix = ".CADD_0.varCount.txt")
prp_r.1k.cadd15 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_1k_Z3_",suffix = ".CADD_15.varCount.txt")
prp_r.1k.cadd25 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_1k_Z3_",suffix = ".CADD_25.varCount.txt")

prp_r.1k.cadd0$DIST <- "1Kb"
prp_r.1k.cadd15$DIST <- "1Kb"
prp_r.1k.cadd25$DIST <- "1Kb"
prp_r.1k.cadd0$Zcutoff <- "3"
prp_r.1k.cadd15$Zcutoff <- "3"
prp_r.1k.cadd25$Zcutoff <- "3"
prp_r.1k.cadd0$CADD <- "CADD_0"
prp_r.1k.cadd15$CADD <- "CADD_15"
prp_r.1k.cadd25$CADD <- "CADD_25"

# Z3,2k:CADD_0/15/25
prp_r.2k.cadd0 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_2k_Z3_",suffix = ".CADD_0.varCount.txt")
prp_r.2k.cadd15 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_2k_Z3_",suffix = ".CADD_15.varCount.txt")
prp_r.2k.cadd25 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_2k_Z3_",suffix = ".CADD_25.varCount.txt")

prp_r.2k.cadd0$DIST <- "2Kb"
prp_r.2k.cadd15$DIST <- "2Kb"
prp_r.2k.cadd25$DIST <- "2Kb"
prp_r.2k.cadd0$Zcutoff <- "3"
prp_r.2k.cadd15$Zcutoff <- "3"
prp_r.2k.cadd25$Zcutoff <- "3"
prp_r.2k.cadd0$CADD <- "CADD_0"
prp_r.2k.cadd15$CADD <- "CADD_15"
prp_r.2k.cadd25$CADD <- "CADD_25"

# Z3,10k:CADD_0/15/25
prp_r.10k.cadd0 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_10k_Z3_",suffix = ".CADD_0.varCount.txt")
prp_r.10k.cadd15 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_10k_Z3_",suffix = ".CADD_15.varCount.txt")
prp_r.10k.cadd25 <- counts2props(IN.DIR,MAFs,TYPEs,prefix = "apa_10k_Z3_",suffix = ".CADD_25.varCount.txt")

prp_r.10k.cadd0$DIST <- "10Kb"
prp_r.10k.cadd15$DIST <- "10Kb"
prp_r.10k.cadd25$DIST <- "10Kb"
prp_r.10k.cadd0$Zcutoff <- "3"
prp_r.10k.cadd15$Zcutoff <- "3"
prp_r.10k.cadd25$Zcutoff <- "3"
prp_r.10k.cadd0$CADD <- "CADD_0"
prp_r.10k.cadd15$CADD <- "CADD_15"
prp_r.10k.cadd25$CADD <- "CADD_25"

# combined
dat.plot <- rbind(prp_r.1k.cadd0,prp_r.1k.cadd15,prp_r.1k.cadd25,
				  prp_r.2k.cadd0,prp_r.2k.cadd15,prp_r.2k.cadd25,
				  prp_r.10k.cadd0,prp_r.10k.cadd15,prp_r.10k.cadd25)

dat.plot$DIST <- factor(dat.plot$DIST,levels = c("1Kb","2Kb","10Kb"))
dat.plot$CADD <- factor(dat.plot$CADD,levels = c("CADD_0","CADD_15","CADD_25"))

saveRDS(dat.plot,file="APA_IPA_RV_enrichment.Z3.add_pvalue.indels.RDS")
#library(ggsci)

#p <- ggplot(dat.plot,aes(x=MAF,y=ESTIM,ymin=CI.LOW,ymax=CI.HIGH,color=DIST)) + 
#  geom_pointrange(linewidth=1,position = position_dodge(width = .7),size=.8) +  
#  geom_hline(yintercept = 1.0, linetype="dashed") + theme_pubr() + scale_color_npg() + 
#  facet_wrap(~Zcutoff,scales = "free_y");p
#p + scale_color_npg()
