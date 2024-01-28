#!/bin/bash

main(){
#	run_define_genebody_region
#	run_extract_genebody_for_outliers
#	run_add_genebody_to_aoutliers
#	run_pick_rv_apa
#	run_liftover
#	run_generate_random_rv
#	run_count_rv_in_RBP_peaks
#	run_enrichment
}

function run_enrichment(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
	for RBP in ${dir}/output/count_rv_in_peaks/*_medz3.real.bed
	do
		RBP_name=`basename $RBP _medz3.real.bed`
		echo $RBP_name
		peak=$(wc -l $RBP | cut -d" " -f 1)
		shuffle_peak=$(wc -l ${dir}/output/count_rv_in_peaks/${RBP_name}_medz3.shuffle.bed |cut -d " " -f 1)
		qtl=$(wc -l ${dir}/output/aOutliers_medz_3.picked_rv.hg19.bed|cut -d" " -f 1)
		echo -e ${RBP_name}"\t"${peak}"\t"${shuffle_peak}"\t"${qtl} >> ${dir}/output/enrichment_genebody_medz3.txt
	done

}

function run_define_genebody_region(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
	cat $dir/input/hg38_Refseq_whole_gene202302328.bed |awk 'BEGIN{OFS="\t"}{if($6=="+" && $3-$2>3000) print $1,$2+3000,$3,$4,$2-$1+3000,$6;else if($6=="-" && $3-$2>3000) print $1,$2,$3-3000,$4,$2-$1+3000,$6}' > ${dir}/input/hg38_Refseq_Genebody.bed
}

function run_extract_genebody_for_outliers(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
#	Rscript $dir/src/extract_genebody_for_outlier.R -r $dir/input/hg38_Refseq_Genebody.bed -g $dir/input/apa_medz.picked.Z_3.rm_global_and_chrX.add_tissueCount.txt -o $dir/input/hg38_Refseq_Genebody.medz_3.apa.bed
	Rscript $dir/src/extract_genebody_for_outlier.R -r $dir/input/hg38_Refseq_Genebody.bed -g $dir/input/apa_medz.picked.Z_2.rm_global_chrX_Y.txt -o $dir/input/hg38_Refseq_Genebody.medz_2.apa.bed
}

function run_add_genebody_to_aoutliers(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
#	python $dir/src/add_genebody_to_outliers.py $dir/input/hg38_Refseq_Genebody.medz_2.apa.v2.bed $dir/input/apa_medz.picked.Z_2.rm_global_chrX_Y.txt|sort -k1,1 -k2,2n > $dir/input/aOutliers_medz_2.txt
	python $dir/src/add_genebody_to_outliers.py $dir/input/hg38_Refseq_Genebody.medz_3.apa.v2.bed $dir/input/apa_medz.picked.Z_3.rm_global_and_chrX.add_tissueCount.txt |sort -k1,1 -k2,2n > $dir/input/aOutliers_medz_3.txt
}

function run_liftover(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
	chain=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/input/hg38ToHg19.over.chain.gz
	liftOver $dir/output/aOutliers_medz_2.picked_rv.uniq_sorted.bed $chain $dir/output/aOutliers_medz_2.picked_rv.hg19.bed $dir/output/medz2.unmapped
	liftOver $dir/output/aOutliers_medz_3.picked_rv.uniq_sorted.bed $chain $dir/output/aOutliers_medz_3.picked_rv.hg19.bed $dir/output/medz3.unmapped
}

function run_generate_random_rv(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
	module load bedtools2/2.29.2

	bedtools shuffle -i $dir/output/aOutliers_medz_3.picked_rv.hg19.bed -g $dir/input/GRCh37_size.txt -chrom -seed 100 -incl ${dir}/input/hg19_Refseq_Genebody.bed -noOverlapping |sort -k1,1 -k2,2n > ${dir}/output/aOutliers_medz_3.shuffle_rv.bed 
	bedtools shuffle -i $dir/output/aOutliers_medz_2.picked_rv.hg19.bed -g $dir/input/GRCh37_size.txt -chrom -seed 100 -incl ${dir}/input/hg19_Refseq_Genebody.bed -noOverlapping |sort -k1,1 -k2,2n > ${dir}/output/aOutliers_medz_2.shuffle_rv.bed 
}

function run_count_rv_in_RBP_peaks(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment
	module load bedtools2/2.29.2
	for RBP in `ls ${dir}/input/CLIP-sig/*.bed`
	do
		RBP_name=`basename $RBP -Sig.bed`
		echo $RBP_name
		bedtools intersect -a $RBP -b $dir/output/aOutliers_medz_2.picked_rv.hg19.bed > ${dir}/output/count_rv_in_peaks/${RBP_name}_medz2.real.bed
		bedtools intersect -a $RBP -b $dir/output/aOutliers_medz_2.shuffle_rv.bed > ${dir}/output/count_rv_in_peaks/${RBP_name}_medz2.shuffle.bed
		bedtools intersect -a $RBP -b $dir/output/aOutliers_medz_3.picked_rv.hg19.bed > ${dir}/output/count_rv_in_peaks/${RBP_name}_medz3.real.bed
		bedtools intersect -a $RBP -b $dir/output/aOutliers_medz_3.shuffle_rv.bed > ${dir}/output/count_rv_in_peaks/${RBP_name}_medz3.shuffle.bed
	done
}

function run_excluding_ov(){

	src=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/src
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/output/APA
	python ${src}/excluding_overlap.py --outliers ${dir}/picked_rv.outliers_allGenes.byGene.bed --controls ${dir}/picked_rv.controls_allGenes.byGene.bed --output ${dir}/picked_rv.outliers_excluding_ov.byGene.bed
	python ${src}/excluding_overlap.py --outliers ${dir}/picked_rv.controls_allGenes.byGene.bed --controls ${dir}/picked_rv.outliers_allGenes.byGene.bed --output ${dir}/picked_rv.controls_excluding_ov.byGene.bed
#	/lustre/home/xdzou/src/anaconda3/bin/python ${src}/excluding_overlap.py --outliers ${dir}/picked_rv.outliers_allGenes.bed --controls ${dir}/picked_rv.controls_allGenes.bed --output ${dir}/picked_rv.outliers_excluding_ov.hg38.bed
#	/lustre/home/xdzou/src/anaconda3/bin/python ${src}/excluding_overlap.py --outliers ${dir}/picked_rv.controls_allGenes.bed --controls ${dir}/picked_rv.outliers_allGenes.bed --output ${dir}/picked_rv.controls_excluding_ov.hg38.bed
}

function run_extract_beta(){
	src=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/src
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects
	python ${src}/integration_with_beta.py --gwas ${dir}/input/UKBB_BMI_bothSex.rareRV_effect.txt --outliers ${dir}/output/GeneExpression/picked_rv.outliers_excluding_ov.hg19.bed --output ${dir}/output/GeneExpression/beta_of_RVs_in_outliers.bed
	python ${src}/integration_with_beta.py --gwas ${dir}/input/UKBB_BMI_bothSex.rareRV_effect.txt --outliers ${dir}/output/GeneExpression/picked_rv.controls_excluding_ov.hg19.bed --output ${dir}/output/GeneExpression/beta_of_RVs_in_controls.bed

	cat ${dir}/output/GeneExpression/beta_of_RVs_in_outliers.bed | awk '$5!="NA"'| sed 's/:/\t/g' > ${dir}/output/GeneExpression/beta_of_RVs_in_outliers.txt &
	wait
	cat ${dir}/output/GeneExpression/beta_of_RVs_in_controls.bed | awk '$5!="NA"'| sed 's/:/\t/g' > ${dir}/output/GeneExpression/beta_of_RVs_in_controls.tmp.bed

	python ${src}/select_controls.py ${dir}/output/GeneExpression/beta_of_RVs_in_outliers.txt ${dir}/output/GeneExpression/beta_of_RVs_in_controls.tmp.bed > ${dir}/output/GeneExpression/beta_of_RVs_in_controls.txt
	rm ${dir}/output/GeneExpression/beta_of_RVs_in_controls.tmp.bed
}

function run_mergeRV_allGene(){
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/output/APA/pick_rv_by_gene
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/output/APA
	chain=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/input/hg38ToHg19.over.chain.gz
	echo "Merge outliers associated RVs ..."
	cat ${indir}/*_outlier.uniq.txt | awk 'BEGIN{OFS="\t"}{fourth=$5":"$6":"$7":"$8":"$1;print $2,$3,$4,fourth}'|sort -k1,1 -k2,2n > ${outdir}/picked_rv.outliers_allGenes.byGene.bed &
	wait
	echo "Merge controls associated RVs ..."
	cat ${indir}/*_control.uniq.txt | awk 'BEGIN{OFS="\t"}{fourth=$5":"$6":"$7":"$8":"$1;print $2,$3,$4,fourth}'|sort -k1,1 -k2,2n > ${outdir}/picked_rv.controls_allGenes.byGene.bed &
	wait

	echo "liftover..."
	liftOver ${outdir}/picked_rv.outliers_allGenes.byGene.bed $chain ${outdir}/picked_rv.outliers_allGenes.hg38Tohg19.byGene.bed ${outdir}/unmapped_outliers.byGene.bed &
	wait
	liftOver ${outdir}/picked_rv.controls_allGenes.byGene.bed $chain ${outdir}/picked_rv.controls_allGenes.hg38Tohg19.byGene.bed ${outdir}/unmapped_controls.byGene.bed &


}


# pick rare variants for each gene in outliers and controls separately
function run_pick_rv_apa(){
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/input/individual_rv
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-03-28-RBP-enrichment/output
	curr_dir=`pwd`

	echo "#!/bin/bash" > ${curr_dir}/submit_pick_rv.slurm
	echo "
#SBATCH --job-name=pickrv
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=pickrv.err
#SBATCH --output=pickrv.out

#################################
module load bedtools2/2.29.2
INDIR=$indir
OUTDIR=$outdir
outlier=$curr_dir/input/aOutliers_medz_3.txt
" >> ${curr_dir}/submit_pick_rv.slurm

	echo $'
cd $SLURM_SUBMIT_DIR

for line in `cat $outliers|awk \'BEGIN{OFS=":"}{print $1,$2,$3,$4,$5,$6}\'`
do
	CHR=`echo $line|awk -F":" \'{print $1}\'`
	S=`echo $line|awk -F":" \'{print $2}\'`
	E=`echo $line|awk -F":" \'{print $3}\'`
	INDS=`echo $line|awk -F":" \'{print $5}\'`
	GENE=`echo $line|awk -F":" \'{print $4}\'`
	echo -e "${CHR}\t${S}\t${E}\t$GENE\t$INDS"| bedtools intersect -a stdin -b ${INDIR}/${INDS}_SNPs_RV.bed | cut -f4- >> ${OUTDIR}/aOutliers_medz_3.picked_rv.txt
done
' >> ${curr_dir}/submit_pick_rv.slurm
#	sbatch ${curr_dir}/submit_pick_rv.slurm 
}


# pick outliers and non-outliers for each gene: one gene one file
function run_pick_aOutliers_and_controls(){
	currDir=`pwd`
	echo "#!/bin/bash" > $currDir/submit_pick_outliers_controls.slurm
	echo "
#SBATCH --job-name=pick_outliers_controls
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=pick_outliers_controls.err
#SBATCH --output=pick_outliers_controls.out

geneBody=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-30-RV-Annotation-VEP/input/All_apa_gene_region.union.sorted.bed
#	outlier=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/output/apa/medz/apa_medz.picked.Z_3.rm_global.txt
outlier=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/input/apa_medz.picked.Z_2.rm_Z3_global_chrXY.txt
medz=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/output/apa/medz/medz_of_all_gene.apa.txt
src=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/src
outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2022-02-17-RV-GWAS-effects/output/APA/pick_outliers_and_controls_by_gene_Z2
" >> $currDir/submit_pick_outliers_controls.slurm
	echo '
cd $SLURM_SUBMIT_DIR

python ${src}/pick_outlier_controls.APA.py --outliers ${outlier} --all_medz $medz --gene_loc ${geneBody} --out_prefix ${outdir}/aoutliers_and_controls_of
' >> $currDir/submit_pick_outliers_controls.slurm

	sbatch $currDir/submit_pick_outliers_controls.slurm
}

main
