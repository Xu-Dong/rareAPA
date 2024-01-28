rm(list = ls())
library(optparse)
library(doParallel)
library(foreach)
library(data.table)
library(dplyr)
library(ggplot2)
library(withr)
library(stringr)
group = "08_5000_500"
if (!file.exists(file.path(group))) {
  dir.create(file.path(group))
}
setwd(group)
num_cores <- 30
registerDoParallel(cores = num_cores)

###############################################################################
files <- list.files(path = "../../2024_01_16_matched_rv/", pattern = "*txt$")
for(m in 1:length(files)){
  flag = TRUE
  tryCatch({
    sample.name <- str_extract(files[m], ".*(?=\\.matched)")
    cat(paste0("Processing ", sample.name, "\n"))
    
    ###############################################################################
    dat <- fread(paste0("../../2024_01_16_matched_rv/", files[m]))
    colnames(dat) <- c("chrom", "pos0", "pos1", "outlier", "aWatershed_score",
                       "Z", "rsChrom", "rsPos0", "rsPos1", "Ref", "Alt", "MAF",
                       "beta", "P", "Quantile")
    # dat <- fread("../2023.08.18_data/1180.matched_rv.rm_indel.tab")
    if(nrow(dat) < 5000){
      next
    } else{
      split_gene <- function(x){
        return(strsplit(x,split=":",fixed=T)[[1]][2])
      }
      
      split_inds <- function(x){
        return(strsplit(x,split=":",fixed=T)[[1]][1])
      }
      
      dat$Gene <- sapply(dat$outlier,split_gene)
      dat$Inds <- sapply(dat$outlier,split_inds)
      # dat <- dat[,c(1:6,11:14,15,16)]
      
      dat$isOutlier[abs(dat$aWatershed_score)>0.8] <- "1"
      dat$isOutlier[abs(dat$aWatershed_score)<0.1] <- "0"
      
      x <- dat %>% filter(isOutlier=="1")
      # x <- x[with_seed(2, sample(1:nrow(x), 40)),]
      outlier_genes <- unique(x$Gene)
      
      dat.f <- dat %>% filter(Gene %in% outlier_genes)
      rm(x,dat)
      dat.f <- na.omit(dat.f)
      dat.f <- distinct(dat.f)
      if(nrow(dat.f) < 500){
        next
      }
      
      ############################################################################
      a <- 0
      b <- 0
      iteration <- 1000
      odds1 <- vector("numeric", length = iteration)
      
      results <- foreach(j = 1:iteration, .packages = c("dplyr")) %dopar% {
        a_local <- 0
        b_local <- 0
        
        cat("permutation: ", j, "\n")
        for (i in 1:length(outlier_genes)) {
          gene <- outlier_genes[i]
          
          dat.g <- dat.f %>% filter(Gene == gene)
          dat.g.o <- dat.g %>% filter(isOutlier == "1")
          dat.g.c <- dat.g %>% filter(isOutlier == "0")
          if (dim(dat.g.o)[1] > 0 && dim(dat.g.c)[1] > 0) {
            set.seed(i)
            dat.g.o <- dat.g.o[sample(1:nrow(dat.g.o), 1), ]
            dat.g.c$ukb_maf_outlier <- dat.g.o$MAF[1]
            dat.g.c$delta_maf <- abs(dat.g.c$MAF - dat.g.c$ukb_maf_outlier)
            dat.g.c <- dat.g.c %>% filter(delta_maf < 1e-3)
            if (dim(dat.g.c)[1] > 0) {
              set.seed(i+j)
              dat.g.c <- dat.g.c[sample(1:nrow(dat.g.c), 1), ]
              beta_o <- abs(dat.g.o$beta)
              beta_c <- abs(dat.g.c$beta)
              
              if (beta_o > beta_c) {
                a_local <- a_local + 1
              } else {
                b_local <- b_local + 1
              }
            }
          }
          if (!(dim(dat.g.o)[1] > 0 && dim(dat.g.c)[1] > 0)) {
            cat(paste0(sample.name, " not meet the standard \n"))
            next
          }
        }
        
        cat("a: ", a_local, " b: ", b_local, "\n")
        if (b_local != 0) {
          return(a_local / b_local)
        } else {
          return(NA)
        }
      }
      odds1 <- unlist(results)
      median(odds1)
      
      ################################################################################
      
      all.matched.nonoutlier <- c()
      for(i in 1:length(outlier_genes)){
        gene.tmp <- outlier_genes[i]
        outlier.tmp <- dat.f %>% filter(isOutlier == 1 & Gene == gene.tmp)
        if(nrow(outlier.tmp) == 0){
          next;
        }
        for(j in 1:nrow(outlier.tmp)){
          outlier.maf <- outlier.tmp[j, ]$MAF
          nonoutlier.tmp <- dat.f %>% filter(isOutlier == 0 
                                             & Gene == gene.tmp 
                                             & MAF >= outlier.maf - 1e-03
                                             & MAF <= outlier.maf + 1e-03)
          
          all.matched.nonoutlier <- unique(append(all.matched.nonoutlier, nonoutlier.tmp$outlier))
        }
      }
      dat.all.matched.nonoutlier <- dat.f[which(dat.f$outlier %in% all.matched.nonoutlier
                                                & dat.f$isOutlier == 0),]
      
      
      #################################################################################
      
      odds2 <- vector("double", length = iteration)
      
      
      results <- foreach(ii = 1:iteration, .packages = c("dplyr")) %dopar% {
        
        aa <- 0
        bb <- 0
        
        
        for (jj in 1:length(unique(dat.all.matched.nonoutlier$Gene))) {
          gene.tmp <- unique(dat.all.matched.nonoutlier$Gene)[jj]
          matched.nonoutlier.gwas.tmp <- dat.all.matched.nonoutlier %>% filter(Gene == gene.tmp)
          
          if (nrow(matched.nonoutlier.gwas.tmp) < 2) {
            next;
          } else {
            variant.select <- matched.nonoutlier.gwas.tmp[with_seed(jj+1, sample(1:nrow(matched.nonoutlier.gwas.tmp), 2)), ]
            
            if(runif(1) < 0.5){
              pseudo.outlier.beta.tmp <- variant.select$beta[1]
              nonoutlier.beta.tmp <- variant.select$beta[2]
            }else{
              pseudo.outlier.beta.tmp <- variant.select$beta[2]
              nonoutlier.beta.tmp <- variant.select$beta[1]
            }
          }
          
          if (abs(pseudo.outlier.beta.tmp) > abs(nonoutlier.beta.tmp)) {
            aa <- aa + 1
          } else {
            bb <- bb + 1
          }
          
        }
        
        aa / bb
      }
      odds2 <- unlist(results)
      median(odds2)
      
      ################################################################
      
      data <- rbind(data.frame(Variant_type = rep("outlier", length(odds1)), Odds = odds1), 
                    data.frame(Variant_type = rep("nonoutlier", length(odds2)), Odds = odds2))
      data <- data[complete.cases(data),]
      test <- wilcox.test(odds1, odds2)
      p.value <- signif(test$p.value, 3)
      p <- ggplot(data, aes(Odds, fill = Variant_type, ..count..))+
        geom_density(alpha = 0.5, adjust = 2.5) +
        labs(y = "Density")+
        theme_bw()+
        theme(plot.title = element_text(hjust = 0.5, size = 12),
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12))+
        scale_fill_manual(values = c("nonoutlier" = "grey", "outlier" = "#3572CB"))+
        scale_x_continuous(limits = c(0,4), breaks = c(0,1,2,3,4))+
        geom_vline(xintercept = c(median(odds2), median(odds1)), color = c("black", "#003f8f"), linetype = "dashed")+
        ggtitle(label = paste0("P-value = ", p.value))
      ggsave(paste0(sample.name, ".density.pdf"), p, width = 6, height = 4.5)
      # print(p)
      
      saveRDS(data, paste0(sample.name, ".odds.rds"))
    }
  }, error = function(e) {
    cat(paste0("Error in ", sample.name, ": ", conditionMessage(e), "\n"))
    flag <- FALSE
  })
  if(!flag) next
}

stopImplicitCluster()
