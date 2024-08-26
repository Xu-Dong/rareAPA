[![RareAPA](https://img.shields.io/badge/RareAPA-v2.1-brightgreen)](https://github.com/Xu-Dong/rareAPA/releases/tag/v2.1)
[![python Release](https://img.shields.io/badge/python-3.8-brightgreen)](https://www.python.org/downloads/)
[![R Release](https://img.shields.io/badge/R-4.3.2-blue)](https://cran.r-project.org/)
![system type](https://img.shields.io/badge/GNU-Linux-brightgreen)
[![RareAPA-Zenodo](https://img.shields.io/badge/Zenodo-blue)](https://doi.org/10.5281/zenodo.10576656)

[Impact of Rare Non-coding Variants on Human Diseases through Alternative Polyadenylation Outliers](https://www.researchsquare.com/article/rs-3907149/v1)

Although rare non-coding variants (RVs) play crucial roles in human complex traits and diseases, understanding their functional mechanisms and identifying those most closely associated with diseases continue to be major challenges. Here, we constructed the first comprehensive atlas of alternative polyadenylation (APA) outliers (aOutliers) from 15,201 samples across 49 human tissues. Strikingly, these aOutliers exhibit unique characteristics markedly distinct from those of outliers based on transcriptional abundance or splicing. This is evidenced by a pronounced enrichment of RVs specifically within aOutliers. Mechanistically, aOutlier RVs frequently alter poly(A) signals and splicing sites, and experimental perturbation of these RVs indeed triggers APA events. Furthermore, we developed a Bayesian-based APA RV prediction model, which successfully pinpointed a specific set of RVs with significantly large effect sizes on complex traits or diseases. A particularly intriguing discovery was the observed convergence effect on APA between rare and common cancer variants, exemplified by the combinatorial regulation of APA in the DDX18 gene. Together, this study introduces a novel APA-enhanced framework for individual genome annotation and underscores the importance of APA in uncovering previously unrecognized functional non-coding RVs linked to human complex traits and diseases.

# Overview

This repository contains codes for data analyses in above rare APA manuscript. Most of the codes are built by R or Python (version Python3) or Shell and all request packages can be found within the codes. Scripts for analyses, including aOutlier calling, rare variants enrichment, RNA Binding Protein (RBP) region enrichment, colocalization with GWAS summary, and aWatershed, can be found in the [Scripts](https://github.com/Xu-Dong/rareAPA/tree/main/scripts) directory.

# Setting up the environment

* R (version > 3.6)
* Python (version 3.8)
* optparse (R package)
* data.table (R package)
* reshape2 (R package)
* plyr/dplyr/magrittr (R packages)
* doMC (R package)
* doParallel (R package)
* foreach (R package)
* withr (R package)
* stringr (R package)
* numpy/pandas (v1.23.3/v1.4.4;python packages)
# Analyses
* APA quantification
  
  APA quantification from multiple samples was conducted by Dapars2. Please find the detailed document through this [linke](https://github.com/3UTR/DaPars2) or this [link](https://github.com/3UTR/3aQTL-pipe)
  
* aOutlier calling

We analyzed single-tissue aOutliers and multi-tissue aOutliers (across at least five tissues). Two scripts, "call_outliers_medz_peerless.xdzou.v3.R" and "call_outliers_single_tissue.py" in the [aOutlier_calling](https://github.com/Xu-Dong/rareAPA/tree/main/scripts/aOutlier_calling) directory perform single-tissue and multi-tissue aOutliers, respectively. Both scripts require one data table containing the profile of normalized APA quantification across individuals and tissues, a demo data can be found in the [Demo]() directory.

* Rare variants enrichment

  Required inputs for RV enrichment analysis includes genotype file of GTEx individuals (obtained from dbGap with request application), variants annotations (VEP and CADD).

* RNA binding region enrichment
  
  Required inputs include RNA binding protein (RBP) binding regions obtained from ENCODE, 166 RBPs were involved, and peak bed file were used.
  
* aOutlier gene enrichemnt with 3'aQTL gene
* Colocalization with GWAS summary data in UK biobank

  
# License

This project is covered under the MIT License.
