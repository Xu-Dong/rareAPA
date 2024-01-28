#!/bin/bash
dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier



chrom=(1 2)
for CHROM in ${chrom[@]}
do
	echo $CHROM
	tabix $dir/input/InDels.tsv.gz $CHROM |cut -f1-4,6 | less > $dir/intput/tmp/Indel_cadd.${CHROM}.tsv

done
