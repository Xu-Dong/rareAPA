#!/usr/bin/env Rscript

####################### FUNCTIONS

#### Function to generate default ggplot colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

#### Function to read in outliers and controls with features for given outlier detection method and variant type
#expression_10k_features_SNPs_MAF0-1.peer_top5.txt
read_type <- function(DIR, MAFs, TYPE, mafstr, prefix){
    pref = paste0(prefix, 'features_')
    ocs = read.table(paste(DIR, pref, TYPE, mafstr, MAFs[1], '.peer_top5.txt', sep = ''), header = T, sep = '\t', stringsAsFactors = F)
    ocs$MAF = MAFs[1]
    if (length(MAFs) > 1) {
        for(i in 2:length(MAFs)){
            ocs_temp = read.table(paste(DIR, pref, TYPE, mafstr, MAFs[i], '.peer_top5.txt', sep = ''),
                                  header = T, sep = '\t', stringsAsFactors = F)
            ocs_temp$MAF = MAFs[i]
            ocs = rbind(ocs, ocs_temp)
        }
    }
    return(ocs)
}

#### Function to read in outliers and controls with features for given outlier detection method across MAFs and variant types
read_file <- function(DIR, MAFs, TYPEs, mafstr, type.names, prefix = ''){
    ocs = read_type(DIR, MAFs, TYPEs[1], mafstr, prefix)
    ocs$TYPE = type.names[1]
    if (length(TYPEs) > 1) {
        for(i in 2:length(TYPEs)){
            ocs_temp = read_type(DIR, MAFs, TYPEs[i], mafstr, prefix)
            ocs_temp$TYPE = type.names[i]
            ocs = rbind(ocs, ocs_temp)
        }
    }
    return(ocs)
}

#### Function to perform logistic regression for each feature against the outlier response
## Input: feature matricx with outlier response
## Output: matrix with odds ratio estimates and p-values
single_logit_core <- function(data, scale = T){
    features = names(data)[!(names(data) %in% c('Y', 'MAF', 'TYPE', 'RANK', 'ID', 'gene_id'))]
    estims = c()
    feature.names = c()
    for(feature in features) {
        if(var(data[, feature], na.rm = TRUE) > 0){
            if(scale){
                lmod = summary(glm(Y ~ scale(data[, feature]), data = data, family = binomial))
            }
            else{
                lmod = summary(glm(Y ~ data[, feature], data = data, family = binomial))
            }
            estims = rbind(estims, matrix(c(lmod$coefficients[2, 1], lmod$coefficients[2, 2], lmod$coefficients[2, 4]), nrow = 1))
            feature.names = append(feature.names, feature)
        }
    }
    estims = data.frame(estims)
    estims = cbind(feature.names, estims)
    names(estims) = c('FEAT', 'BETA', 'STD', 'PVAL')
    return(estims)
}

#### Function to calculate functional enrichments for given outlier detection method and variant type across MAFs
single_logit_type <- function(data, MAFs, TYPE, scale = T){
    estims = single_logit_core(data[data$MAF == MAFs[1] & data$TYPE == TYPE, ], scale = scale)
    estims$MAF = MAFs[1]
    estims$TYPE = TYPE
    if (length(MAFs) > 1) {
        for(i in 2:length(MAFs)){
            estims_temp = single_logit_core(data[data$MAF == MAFs[i] & data$TYPE == TYPE, ], scale = scale)
            estims_temp$MAF = MAFs[i]
            estims_temp$TYPE = TYPE
            estims = rbind(estims, estims_temp)
        }
    }
    return(estims)
}

#### Function to calculate functional enrichments for given outlier detection method across variant types and MAFs
single_logit <- function(data, MAFs, TYPEs, scale = T){
    estims = single_logit_type(data, MAFs, TYPEs[1], scale = scale)
    if (length(TYPEs) > 1) {
        for(i in 2:length(TYPEs)){
            estims_temp = single_logit_type(data, MAFs, TYPEs[i], scale = scale)
            estims = rbind(estims, estims_temp)
        }
    }
    return(estims)
}

#### Function for scientific notation 
fancy_scientific <- function(l) {
    ## turn in to character string in scientific notation
    lnew <- format(l, scientific = TRUE)
    ## quote the part before the exponent to keep all the digits
    lnew <- gsub("^(.*)e", "'\\1'e", lnew)
    ## turn the 'e+' into plotmath format
    lnew <- gsub("e", "%*%10^", lnew)
    ## don't use scientific notation if exponent is zero
    for (i in 1:length(lnew)) {
        if (length(grep('+00', lnew[i], fixed = TRUE)) > 0) {
            lnew[i] = l[i]
        }
    }
    ## return this as an expression
    return(parse(text=lnew))
}


#### for each variant type and each minor allele frequency bin,
#### get the ratio of the portions of outliers/non outliers with rare variants
#### and the associated confidence intervals (exp of Wald interval on the log of the proportion ratio)
proportion.ratios = function(counts) {
    mafs = unique(counts$MAF)
    types = unique(counts$TYPE)
    ## make emty data frame and fill it from rbinds (not very many, so it's fine)
    results = data.frame(ESTIM=numeric(), CI.LOW=numeric(), CI.HIGH=numeric(), MAF=character(), TYPE=character(),
                         stringsAsFactors = F)
    for (m in mafs) {
        for (t in types) {
            results = rbind(results, proportion.ratios.helper(counts, m, t))
        }
    }
    return(results)
}

#### actually does the work of the main function for the given maf bin and variant type
proportion.ratios.helper = function(counts, maf, type) {
    counts.subset = counts[counts$MAF == maf & counts$TYPE ==type, c('Y','any_variant')]
    summary.counts = as.data.frame(table(counts.subset))
    stopifnot(nrow(summary.counts)== 4 & min(summary.counts$Freq) > 0) # check that no zero values, and right number of values
    ## get required values to get statistic and CI
    out.var = summary.counts$Freq[summary.counts$Y == 1 & summary.counts$any_variant == 1]
    nonout.var = summary.counts$Freq[summary.counts$Y == 0 & summary.counts$any_variant == 1]
    out.total = sum(summary.counts$Freq[summary.counts$Y == 1])
    nonout.total = sum(summary.counts$Freq[summary.counts$Y == 0])
    ## actually calculate
    estimate = (out.var/out.total)/(nonout.var/nonout.total)
    ## get bounds of confidence interval on the log of the proportion then exponentiate
    log.se = sqrt(1/out.var - 1/out.total + 1/nonout.var - 1/nonout.total)
    max.ci = estimate * exp(1.96*log.se)
    min.ci = estimate * exp(-1.96*log.se)
    ## put all together in a list that can become the row of a dataframe
    dfrow = list(ESTIM=estimate, CI.LOW=min.ci, CI.HIGH=max.ci, MAF=maf, TYPE=type)
    return (dfrow)
}
