#!/bin/bash

main(){
#	count_RV_by_anno 1k 1000 PAS50bp SNP
#	count_RV_by_anno 1k 1000 PAS50bp indel
#	count_RV_by_anno 1k 1000 Coding SNP
#	count_RV_by_anno 1k 1000 Conserved SNP
#	count_RV_by_anno 1k 1000 UTR_3 SNP
#	count_RV_by_anno 1k 1000 Splice SNP

#	run_features Coding SNP

#	run_features Conserved SNP
#	sleep 1800
	run_features Splice SNP
	sleep 1800
	run_features UTR_3 SNP
	sleep 1800
	run_features PAS50bp SNP

#	run_summarize_and_enrichment
}


# functions

#function preprocess_indel_cadd(){
#}

function run_merge_snp_indel_count(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
#	ANNO=(Coding Conserved Enhancer Frame_shift Intergenic Intron Promoter Splice Stop TFBS UTR_3 UTR5)
	ANNO=(PAS50bp)
	for anno in ${ANNO[@]}
	do
		echo $anno
		outdir=$dir/output/Count_by_Anno/1k/$anno/SNP_indel
		if [ ! -d $outdir ]
		then
			mkdir -p $outdir
		fi

		echo "#!/bin/bash" > $dir/submit_merge_count_${anno}.slurm
		echo "
#SBATCH --job-name=$anno
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=${anno}.err
#SBATCH --output=${anno}.out

DIR=$dir
FEATURE=$anno
###############################
" >> $dir/submit_merge_count_${anno}.slurm
		echo '
echo "process start at:"
date
cd $SLURM_SUBMIT_DIR

Rscript $DIR/src/merge_SNP_indel_count.R -a $FEATURE

' >> $dir/submit_merge_count_${anno}.slurm
		sbatch $dir/submit_merge_count_${anno}.slurm
	done
}

function prepare_RV_annotation_file(){
	varType=$1
	currdir=`pwd`

	for f in `ls $currdir/Task/task_*`
	do
		task=`basename $f`
		echo $task
		echo "#!/bin/bash" > $currdir/submit_prepare_RVanno_${task}.slurm
		echo "
#SBATCH --job-name=Indel_$task
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=indel_${task}.err
#SBATCH --output=indel_${task}.out

##############################
VAR=$varType
CURRDIR=$currdir
TASK=$task
" >> $currdir/submit_prepare_RVanno_${task}.slurm

		echo '
echo "process start at:"
date
cd $SLURM_SUBMIT_DIR

for inds in `cat $CURRDIR/Task/$TASK`
do
	echo $inds
	Rscript $CURRDIR/src/extract_gtex_variant_PASanno.R -i $inds -v $VAR
done

echo "process end at:"
date
' >> $currdir/submit_prepare_RVanno_${task}.slurm
		sbatch $currdir/submit_prepare_RVanno_${task}.slurm
	done
}


function prepare_outliers(){
	dir=`pwd`
	Rscript $dir/src/prepare_outliers.R
}

function count_RV_by_anno(){
	window=$1
	window_size=$2
	anno=$3
	varType=$4
	UTR=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/All_singleTissue_aOutlier_genebody.bed
#	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/GTEx_Variants_with_PASanno
	indir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/GTEx_Variants_with_Anno
	outdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/Count_by_Anno/${window}/${anno}
	current_dir=`pwd`
	if [ ! -d ${outdir} ]; then
		mkdir -p ${outdir}
		cd ${outdir}
		mkdir SNP
		mkdir indel
		cd ..
	fi
	cd $current_dir
	
	echo "#!/bin/bash" > ${current_dir}/submit_countRV_${window}_${anno}.${varType}.slurm
	echo "
#SBATCH --job-name=count_${anno}_${varType}_$window
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error=count_${anno}_${varType}_${window}.err
#SBATCH --output=count_${anno}_${varType}_${window}.out

##################################
module load bedtools2/2.29.2
indir=$indir
outdir=$outdir
UTR=$UTR
window=$window
window_size=$window_size
nproc=32
annotation=$anno
VARtype=$varType
" >> ${current_dir}/submit_countRV_${window}_${anno}.${varType}.slurm

	echo $'
makecountfeatures(){
	var_bed=$1
	fname=`basename $var_bed`
	sample=${fname%%_*}
	vartype=$VARtype
	
	out_prefix=${sample}_counts_bygene
	outfile=${outdir}/${vartype}/${out_prefix}.txt
	echo -e "gene_id\tn_variants" > $outfile
	cat $var_bed | cut -f1-4,13 |\
		awk -v ANNO=$annotation \'$5==ANNO\' | \
	bedtools window -c -r $window_size -l $window_size -a $UTR -b stdin | cut -f4,5 >> $outfile
}

export -f makecountfeatures
export window_size
export UTR
export outdir
export VARtype
export indir
export annotation

echo "Processing SNPs and indels..."
parallel --jobs $nproc makecountfeatures ::: ${indir}/GTEX-*_SNP_anno.bed
#parallel --jobs $nproc makecountfeatures ::: ${indir}/GTEX-*_indel_anno.bed

module unload bedtool2/2.29.2' >> ${current_dir}/submit_countRV_${window}_${anno}.${varType}.slurm
	sbatch ${current_dir}/submit_countRV_${window}_${anno}.${varType}.slurm
}




run_features(){
	anno=$1
	varType=$2
	window=1k
	threshold=5
	maf="0-1"
	
	inputdir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/Count_by_Anno/${window}/$anno/SNP
	inds=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/input/gtex_inds_715.txt
	stattype=Z
	suffix="_counts_bygene.txt"
	curr_dir=`pwd`
	cd $curr_dir

	for tissue in `cat $curr_dir/gtex_49tissues.txt`
	do
		echo $tissue
		outliers=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/Z${threshold}/${tissue}.outlier.txt
		counts=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/Zmatrix/${tissue}.z_mat.txt
		output_file=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/single_tissue/enrichment/varCount_${window}_Z${threshold}_${varType}_MAF${maf}.${anno}.${tissue}.txt

		echo "#!/bin/bash" > $curr_dir/submit_varcount_${tissue}.slurm
		echo "
#SBATCH --job-name=varCount_${tissue}
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=varCount_${tissue}.err
#SBATCH --output=varCount_${tissue}.out

##################################
script=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier/src
inputdir=${inputdir}
inds=$inds
outliers=${outliers}
counts=$counts
threshold=$threshold
stattype=$stattype
output_file=$output_file
suffix=$suffix" >> $curr_dir/submit_varcount_${tissue}.slurm

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
' >> $curr_dir/submit_varcount_${tissue}.slurm
		sbatch $curr_dir/submit_varcount_${tissue}.slurm
	done
}


function run_summarize_and_enrichment(){
	dir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-06-03-RV-enrichment-to-outlier
	Rscript $dir/src/enrichment_RV_fisher.SNV.R
}

# ---------------------------------------------- Main body ---------------------------------------------------------------------------

main
