#!/bin/bash

vcf=/lustre/home/llei/Data/GTEx_v8/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz

# individuals to include from the allele frequency calculation
# generate file name from the indincl file name
keeplist=/lustre/home/xdzou/RareVar_aQTL_Project/GTExV8/individuals_VCFids_EA.txt

# SNPs
#vcftools --gzvcf $vcf --out singleton_SNPs --remove-filtered-all --keep $keeplist --remove-indels --max-missing-count 10 --singletons 

# indels
vcftools --gzvcf $vcf --out singleton_indels --remove-filtered-all --keep $keeplist --keep-only-indels --max-missing-count 10 --singletons

