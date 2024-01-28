#!/bin/bash

main(){
#	get_single_tissue_aoutlier_gene_list
#	add_genebody_region
#	split_aoutlier_by_tissue
#	count_apafeature_utr1k 1k 1000 0
#	count_apafeature_utr1k 1k 1000 15
#	count_apafeature_utr1k 1k 1000 25
#	run_features SNPs 0-1 1k 0
#	run_features SNPs 0-1 1k 15
#	run_features SNPs 0-1 1k 25
#	run_features "SNPs" "-S" "1k" "0"
#	run_features "SNPs" "-S" "1k" "15"
#	run_features "SNPs" "-S" "1k" "25"
#	run_features indels 0-1 1k 0
#	run_features indels 0-1 1k 15
	run_combine_snv_indel
}


# --- functions
function get_single_tissue_aoutlier_gene_list(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier

	cat $dir/input/aOutliers_singlez_picked.Z3.txt|cut -f1|tail -n+2|sort|uniq > $dir/input/single_tissue_aOutlier.Z3.gene_list.txt
}


function add_genebody_region(){
# add genebody region in "All_apa_gene_region.union.sorted.bed" to single-tissue aOutlier gene list
# genes in chrX and chrY were removed from the final list
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
	Rscript $dir/src/add_gb_for_singletissue_aoutlier.R -r $dir/input/All_apa_gene_region.union.sorted.bed -s $dir/input/single_tissue_aOutlier.Z3.gene_list.txt -o $dir/input/All_singleTissue_aOutlier_genebody.bed

}

function split_aoutlier_by_tissue(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
	outlier_dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/output/apa

	Rscript $dir/src/split_singleTissue_outlier_by_tissue.R -i $outlier_dir/aOutliers_singlez_picked.Z6.txt -t $dir/single_tissue/gtex_49tissues.txt

}

function count_apafeature_utr1k(){
	window=$1
	window_size=$2
	cadd=$3
	UTR=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/All_singleTissue_aOutlier_genebody.bed
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/cadd_singleton
#	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/cadd_singleton_indel
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/anno_byGene/${window}/CADD_${cadd}
	current_dir=`pwd`
	if [ ! -d ${outdir}/counts ]; then
		mkdir -p ${outdir}/counts
		cd ${outdir}/counts
		mkdir MAF-S
		mkdir MAF0-1
		for dir in MAF*
		do
			cd $dir
			mkdir SNPs
			mkdir indels
			cd ..
		done
	fi
	cd $current_dir
	
	endpoints=(0 0.01)
	names=(0 1)
	for ((i=1; i<${#endpoints[@]}; i++))
	do
		maf=MAF${names[$i-1]}-${names[$i]}
		echo "#!/bin/bash" > ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.SNP.slurm
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
end_point=${endpoints[$i]}" >> ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.SNP.slurm

		echo $'
makecountfeatures(){
	var_bed=$1
	fname=`basename $var_bed`
	sample=${fname%%_*}
	vartype=indels
	
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

module unload bedtool2/2.29.2' >> ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.SNP.slurm
		sbatch ${current_dir}/submit_countFeatures_${window}_CADD${cadd}_${maf}.SNP.slurm
	done
}




run_features(){
	vartype=$1
	maf=$2
	window=$3
	cadd=$4
	threshold=3
	
	inputdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/anno_byGene/${window}/CADD_${cadd}/counts/MAF${maf}/${vartype}
	inds=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/gtex_inds_715.txt
	stattype=Z
	suffix="_counts_bygene.txt"
	curr_dir=`pwd`
	cd $curr_dir

	for tissue in `cat $curr_dir/gtex_49tissues.txt`
	do
		echo "$tissue"
		outliers=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/Z${threshold}/${tissue}.outlier.txt
		counts=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/Zmatrix/${tissue}.z_mat.txt
		output_file=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/enrichment/varCount_${window}_Z${threshold}_${vartype}_MAF${maf}.CADD_${cadd}.${tissue}.txt

		echo "#!/bin/bash" > $curr_dir/submit_${tissue}_${window}_${maf}_CADD${cadd}.slurm
		echo "
#SBATCH --job-name=${tissue}_$maf_$cadd
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=${tissue}_${maf}_${cadd}.err
#SBATCH --output=${tissue}_${maf}_${cadd}.out

##################################
script=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/src
inputdir=${inputdir}
inds=$inds
outliers=${outliers}
counts=$counts
threshold=$threshold
stattype=$stattype
output_file=$output_file
suffix=$suffix" >> $curr_dir/submit_${tissue}_${window}_${maf}_CADD${cadd}.slurm

		echo '
echo "process start at:"
date

cd $SLURM_SUBMIT_DIR

/lustre/home/xdzou/src/anaconda2/bin/python ${script}/pick_outliers_controls_imbalanced.singleTissue.py --FEATURE_DIR $inputdir \
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
' >> $curr_dir/submit_${tissue}_${window}_${maf}_CADD${cadd}.slurm
		sbatch $curr_dir/submit_${tissue}_${window}_${maf}_CADD${cadd}.slurm
	done
}


function run_summarize_and_enrichment(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
	Rscript $dir/src/enrichment_RV_fisher.SNV.R
}


function run_combine_snv_indel(){
	curr_dir=`pwd`
	
	for tissue in `cat $curr_dir/gtex_49tissues.txt`
	do
		echo "$tissue"
		echo "#!/bin/bash" > $curr_dir/submit_combine_${tissue}.slurm
		echo "
#SBATCH --job-name=${tissue}
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=${tissue}.err
#SBATCH --output=${tissue}.out

##################################
DIR=$curr_dir
TISSUE=$tissue
" >> $curr_dir/submit_combine_${tissue}.slurm

		echo '
echo "process start at:"
date

cd $SLURM_SUBMIT_DIR

Rscript $DIR/merge_snv_indel.R -t $TISSUE


echo "process end at:"
date
' >> $curr_dir/submit_combine_${tissue}.slurm
#		sbatch $curr_dir/submit_combine_${tissue}.slurm
	done

}


# --- 
main
