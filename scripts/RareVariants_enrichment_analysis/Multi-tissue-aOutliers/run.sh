#!/bin/bash
main(){
#	extract_singleton_in_GTEx_EA
#	preprocess_indel_cadd
#	trim_anno_bySite SNPs
#	trim_anno_bySite indels
#	anno_cadd_indel
#	annotate_singleton
#	prepare_outliers
	count_apafeature_utr1k 1k 1000 0
#	count_apafeature_utr1k 1k 1000 15
#	count_apafeature_utr1k 1k 1000 25

#	count_apafeature_utr1k 2k 2000 0
#	count_apafeature_utr1k 2k 2000 15
#	count_apafeature_utr1k 2k 2000 25

#	count_apafeature_utr1k 10k 10000 0
#	count_apafeature_utr1k 10k 10000 15
#	count_apafeature_utr1k 10k 10000 25

#	run_features SNPs 0-1 10k 0
#	run_features SNPs 1-5 10k 0
#	run_features SNPs 5-10 10k 0
#	run_features SNPs 10-25 10k 0
##	run_features "SNPs" "-S" "10k" "0"
#	run_features SNPs 0-1 10k 15
#	run_features SNPs 1-5 10k 15
#	run_features SNPs 5-10 10k 15
#	run_features SNPs 10-25 10k 15
#	run_features "SNPs" "-S" "10k" "15"

#	run_features indels 0-1 10k 15
#	run_features indels 1-5 10k 15
#	run_features indels 5-10 10k 15
#	run_features indels 10-25 10k 15
#	run_features "indels" "-S" "10k" "15"
#	run_features indels 0-1 10k 0
#	run_features indels 1-5 10k 0
#	run_features indels 5-10 10k 0
#	run_features indels 10-25 10k 0
#	run_features "indels" "-S" "10k" "0"
# summarize Z3:1k/2k/10k:CADD_0/15/25
#	run_summarize_and_enrichment
}


# functions

function extract_singleton_in_GTEx_EA(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
	bash $dir/src/get_singletons.sh
}

function preprocess_indel_cadd(){
	currDir=`pwd`

	echo "#!/bin/bash" > $currDir/submit_preprocess_indel_cadd.sh
	echo "
dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
	" >> $currDir/submit_preprocess_indel_cadd.sh
	echo '
for CHROM in ${chrom[@]}
do
	echo $CHROM
	tabix $dir/input/InDels.tsv.gz $CHROM |cut -f1-4,6 | less > $dir/intput/tmp/Indel_cadd.${CHROM}.tsv
done
	'
}

function trim_anno_bySite(){
	varType=$1
	currdir=`pwd`
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-05-05-outlier-gene-annotation/output/anno_bySite
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/GTEx_Variants_by_INDS

	for f in `ls $currdir/Task/task_*|tail -n+3`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash" > $currdir/submit_trim_anno_bySite_${task}.slurm
		echo "
#SBATCH --job-name=$task
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=${task}.err
#SBATCH --output=${task}.out

##############################
VAR=$varType
INDIR=$indir
OUTDIR=$outdir
CURRDIR=$currdir
TASK=$task
CADD=$currdir/input/GTEx_all_indels.cadd.bed
" >> $currdir/submit_trim_anno_bySite_${task}.slurm

		echo '
echo "process start at:"
date
cd $SLURM_SUBMIT_DIR

for inds in `cat $CURRDIR/Task/$TASK`
do
	echo $inds
#	zcat $INDIR/${inds}_${VAR}_features.bed.gz|cut -f1-4,7|sort|uniq > $OUTDIR/${inds}_SNPs_cadd.bed
	zcat $INDIR/${inds}_${VAR}_features.bed.gz|cut -f1-4|sort|uniq |bedtools intersect -wo -a stdin -b $CADD|cut -f1-4,8> $OUTDIR/${inds}_indels_cadd.bed
done

echo "process end at:"
date
' >> $currdir/submit_trim_anno_bySite_${task}.slurm
		sbatch $currdir/submit_trim_anno_bySite_${task}.slurm
	done
}

function anno_cadd_indel(){
	currdir=`pwd`
	chrom_list=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
	for chrom in ${chrom_list[@]}
	do
		echo chr$chrom
		echo "#!/bin/bash" > $currdir/submit_annoCADD_${chrom}.slurm
		echo "
#SBATCH --job-name=$chrom
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=chr_${chrom}.err
#SBATCH --output=chr_${chrom}.out
##############################

dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
chromosome=$chrom
" >> $currdir/submit_annoCADD_${chrom}.slurm
		echo '
cd $SLURM_SUBMIT_DIR
echo "process start at:"
date
Rscript $dir/src/anno_CADD_for_Indels.R -c $chromosome

echo "process end at:"
date
' >> $currdir/submit_annoCADD_${chrom}.slurm
		sbatch $currdir/submit_annoCADD_${chrom}.slurm

	done

#	cat GTEx_indel_cadd.*.txt|grep -v "^chrom"|awk '{print $1"\t"($2-1)"\t"$2"\t"$3}'|sort -k1,1 -k2,2n > ../GTEx_all_indels.cadd.bed
}

function annotate_singleton(){
	currdir=`pwd`
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/GTEx_Variants_by_INDS
#	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/cadd_singleton
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/cadd_singleton_indel
	if [ ! -d "$outdir" ]
	then
		mkdir -p $outdir
	fi

	for f in `ls $currdir/Task/task_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash" > $currdir/submit_anno_singleton_${task}.slurm
		echo "
#SBATCH --job-name=$task
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=${task}.err
#SBATCH --output=${task}.out

##############################
INDIR=$indir
OUTDIR=$outdir
CURRDIR=$currdir
TASK=$task
" >> $currdir/submit_anno_singleton_${task}.slurm

		echo '
echo "process start at:"
date
cd $SLURM_SUBMIT_DIR

for inds in `cat $CURRDIR/Task/$TASK`
do
	echo $inds
#	Rscript $CURRDIR/src/label_singleton_doubleton.R -s $CURRDIR/input/singleton_SNPs.singletons -i $INDIR/${inds}_SNPs_cadd.bed -o $OUTDIR/${inds}_cadd_singleton.bed
	Rscript $CURRDIR/src/label_singleton_doubleton.R -s $CURRDIR/input/singleton_indels.singletons -i $INDIR/${inds}_indels_cadd.bed -o $OUTDIR/${inds}_cadd_singleton.bed
done

echo "process end at:"
date
' >> $currdir/submit_anno_singleton_${task}.slurm
		sbatch $currdir/submit_anno_singleton_${task}.slurm
	done
}

function prepare_outliers(){
	dir=`pwd`
	Rscript $dir/src/prepare_outliers.R
}

function count_apafeature_utr1k(){
	window=$1
	window_size=$2
	cadd=$3
	UTR=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/All_apa_ipa.genebody.bed
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/cadd_singleton
	#indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/cadd_singleton_indel
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/output/anno_byGene/${window}/CADD_${cadd}
	current_dir=`pwd`
	if [ ! -d ${outdir}/counts ]; then
		mkdir -p ${outdir}/counts
		cd ${outdir}/counts
		mkdir MAF-S
#		mkdir MAF-D
		mkdir MAF0-1
		mkdir MAF1-5
		mkdir MAF5-10
		mkdir MAF10-25
		for dir in MAF*
		do
			cd $dir
			mkdir SNPs
			mkdir indels
			cd ..
		done
	fi
	cd $current_dir
	
	endpoints=(0 0.01 0.05 0.1 0.25)
	names=(0 1 5 10 25)
	for ((i=1; i<${#endpoints[@]}; i++))
	do
		maf=MAF${names[$i-1]}-${names[$i]}
		echo "#!/bin/bash" > ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.indel.slurm
		echo "
#SBATCH --job-name=countFeature_CADD${cadd}_${maf}_$window
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error=countFeature_CADD${cadd}_${maf}_${window}.err
#SBATCH --output=countFeature_CADD${cadd}_${maf}_${window}.out

##################################
module load bedtools2/2.29.2
indir=$indir
outdir=$outdir
UTR=$UTR
window=$window
window_size=$window_size
nproc=32
maf=$maf
CADD=$cadd
start_point=${endpoints[$i-1]}
end_point=${endpoints[$i]}" >> ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.indel.slurm

		echo $'
makecountfeatures(){
	var_bed=$1
	fname=`basename $var_bed`
	sample=${fname%%_*}
	vartype=SNPs
	
	out_prefix=${sample}_counts_bygene
	mafname=${maf}
	outfile=${outdir}/counts/${mafname}/${vartype}/${out_prefix}.txt
	echo -e "gene_id\tn_variants" > $outfile
	cat $var_bed | \
		awk -v start=${start_point} -v end=${end_point} -v cutCADD=$CADD \'$5>cutCADD && ($6!="S" || $6!="D") && $4>start && $4<=end\' | \
	bedtools window -c -r $window_size -l $window_size -a $UTR -b stdin | cut -f4,5 >> $outfile
}

export -f makecountfeatures
export window_size
export UTR
export outdir
export CADD
export indir
export maf
export start_point
export end_point

echo "Processing SNPs and indels..."
parallel --jobs $nproc makecountfeatures ::: ${indir}/GTEX-*

module unload bedtool2/2.29.2' >> ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.indel.slurm
#		sbatch ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.indel.slurm
	done
}




run_features(){
	vartype=$1
	maf=$2
	window=$3
	cadd=$4
	threshold=2
	
	inputdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/output/anno_byGene/${window}/CADD_${cadd}/counts/MAF${maf}/${vartype}
	inds=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/gtex_inds_715.txt
	outliers=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/outliers/apa_medz.picked.Z_${threshold}.txt
	counts=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/outliers/apa_medz.counts.txt
	stattype=Z
	output_file=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/output/enrichment/apa_${window}_Z${threshold}_${vartype}_MAF${maf}.CADD_${cadd}.varCount.txt
	suffix="_counts_bygene.txt"
	curr_dir=`pwd`
	cd $curr_dir
	echo "#!/bin/bash" > $curr_dir/submit_varcount_${window}_${maf}_CADD${cadd}.slurm
	echo "
#SBATCH --job-name=varCount_${window}_$maf_$cadd
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=varCount_${window}_${maf}_${cadd}.err
#SBATCH --output=varCount_${window}_${maf}_${cadd}.out

##################################
script=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/src
inputdir=${inputdir}
inds=$inds
outliers=${outliers}
counts=$counts
threshold=$threshold
stattype=$stattype
output_file=$output_file
suffix=$suffix" >> $curr_dir/submit_varcount_${window}_${maf}_CADD${cadd}.slurm

	echo '
echo "process start at:"
date

cd $SLURM_SUBMIT_DIR

/lustre/home/xdzou/src/anaconda2/bin/python ${script}/pick_outliers_controls_imbalanced.py --FEATURE_DIR $inputdir \
--OUTLIER_PICKED $outliers \
--COUNTS $counts \
--INDS $inds \
--OUT $output_file \
--type $stattype \
--threshold $threshold \
--END $suffix &
wait

echo "process end at:"
date
' >> $curr_dir/submit_varcount_${window}_${maf}_CADD${cadd}.slurm
	sbatch $curr_dir/submit_varcount_${window}_${maf}_CADD${cadd}.slurm
}


function run_summarize_and_enrichment(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
	Rscript $dir/src/enrichment_RV_fisher.SNV.R
}

# ---------------------------------------------- Main body ---------------------------------------------------------------------------

main
