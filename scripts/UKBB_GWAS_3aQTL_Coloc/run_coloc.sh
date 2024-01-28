#!/bin/bash

# --------------------- Functions ----------------------
# coloc_5.1.0.1
# main function

main(){
	#prepare_trait_list
	#run_coloc /lustre/home/lhgong/2023-04-21-UKB-coloc/input/trait_list.txt
	#find_colocalized_gene_all
	#find_sigSentinalSnp_all
}


prepare_trait_list() {
	working_dir=`pwd`
	GWAS_dir="/lustre/home/lhgong/share/ukbb_coloc/coloc_format_data"
	
	ls /lustre/home/lhgong/share/ukbb_coloc/coloc_format_data | cut -d . -f1 > ${working_dir}/input/trait_list.txt

	mkdir -p ${working_dir}/input/gwas_dir
	for gwas in `cat ${working_dir}/input/trait_list.txt`
	do
		gwas_file="/lustre/home/lhgong/share/ukbb_coloc/coloc_format_data/${gwas}.coloc_format.txt"
		mkdir -p ${working_dir}/input/gwas_dir/${gwas}
		echo "ln -s $gwas_file ${working_dir}/input/gwas_dir/${gwas}/GLGC_${gwas}_results.txt"
		ln -s $gwas_file ${working_dir}/input/gwas_dir/${gwas}/GLGC_${gwas}_results.txt
	done
}

function find_colocalized_gene_all(){
	dir=/lustre/home/lhgong/2023-04-21-UKB-coloc
	if [ ! -d "$dir/output/summary_coloc" ]
	then
		mkdir -p $dir/output/summary_coloc
	fi
	python $dir/src/find_colocalizedGene_v2_li.py
}

function find_sigSentinalSnp_all(){
	dir=/lustre/home/lhgong/2023-04-21-UKB-coloc

	if [ ! -d "$dir/output/summary_coloc" ]
	then
		mkdir -p $dir/output/summary_coloc
	fi

	python $dir/src/sigSentinalSnp_all_li.py
}

function run_coloc(){
# coloc version: coloc_5.1.0.1
	var=1
	task=$1
	#GWASs_hg38=/lustre/home/lhgong/2023-04-21-UKB-coloc/input/UKB_GWAS_coloc_format
	GWAS_dir="/lustre/home/lhgong/2023-04-21-UKB-coloc/input/gwas_dir"
	QTLfolder=/lustre/home/hchen/2021-10-31-cancer-GWAS/aQTL_coloc/output/aQTLs_transcript_hg19
	srcFolder=/lustre/home/lhgong/2023-04-21-UKB-coloc/src/src_pvalue_maf

	curr_dir="/lustre/home/lhgong/2023-04-21-UKB-coloc"
	echo ${task}

	for trait in `cat $task`
	do
		echo $var
		echo $GWAS_dir/${trait}/GLGC_${trait}_results.txt

		if [ -f $GWAS_dir/${trait}/GLGC_${trait}_results.txt ]
		then
			echo -e "$task : $trait"
			for folder in $QTLfolder/*
			do
				if [ "`ls $folder`" != "" ]
				then
					tissue=`echo "$folder" | awk -F"/" '{print $NF }'`
					mkdir -p ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}
					if [ -d "${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/src_pvalue_maf/" ]
					then
						rm -rf ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/src_pvalue_maf/ &
						wait
					fi
					cp -r $srcFolder/ ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait} &
					wait
					cd ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}

					if [[ -f "summary_table.txt" ]] && [[ -f "coloc_bed_table.BED" ]]
					then
						echo -e "$tissue\t$trait exists"
						continue
					else
						cd ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/src_pvalue_maf/
						echo "#!/bin/bash" > ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/submit_coloc.slurm
						echo "
#SBATCH --job-name=${trait}_$tissue
#SBATCH --partition=cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/${trait}_${tissue}.err
#SBATCH --output=${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/${trait}_${tissue}.out

DIR=${curr_dir}
TISSUE=$tissue
TRAIT=$trait
GWAS_path=$GWAS_dir
QTL_path=$QTLfolder
Rpath=/lustre/home/wychen/anaconda3/envs/r4/bin" >>${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/submit_coloc.slurm

						echo $'
echo "process will start at:"
date
echo "+++++++++++++++++++++++++++"
module load intel_mpi/2018_u1
module load R/3.6.2-anaconda3

cd ${DIR}/output/coloc_res/${TRAIT}/${TISSUE}_${TRAIT}/src_pvalue_maf

${Rpath}/Rscript ${DIR}/output/coloc_res/${TRAIT}/${TISSUE}_${TRAIT}/src_pvalue_maf/wrapper.R ./ ${QTL_path}/${TISSUE}/ $GWAS_path/${TRAIT}/ /lustre/home/lhgong/share/ukbb_coloc/sentinal_data/${TRAIT}.sentinal.txt &
wait

mv summary_table.txt coloc_bed_table.BED ../
wait
echo "process will sleep 30s"
sleep 30
echo "process will end at:"
date

module unload intel_mpi/2018_u1
module unload R/3.6.2-anaconda3
' >> ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/submit_coloc.slurm &
						#wait
						#sbatch ${curr_dir}/output/coloc_res/${trait}/${tissue}_${trait}/submit_coloc.slurm &
						wait
					fi
				fi

			done
			echo "$trait finished!"
		fi
	done
}

function run_format_qtl(){
	qtl_dir=/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project/output/QTL_mapping
	curr_dir=`pwd`
	if [ ! -d "$curr_dir/input/cis_QTL" ]
	then
		mkdir -p $curr_dir/input/cis_QTL
	fi

	while read line
	do
		tissue=`echo $line|awk '{print $1}'`
		size=`echo $line|awk '{print $2}'`
		echo "Current tissue and size: $tissue,$size"
		echo "#!/bin/bash" > ${curr_dir}/submit_formatQTL_${tissue}.slurm
		echo "
#SBATCH --job-name=format_$tissue
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=format_${tissue}.err
#SBATCH --output=format_${tissue}.out

QTLdir=$qtl_dir
TISSUE=$tissue
N=$size
baseDir=$curr_dir
" >> ${curr_dir}/submit_formatQTL_${tissue}.slurm
		echo $'
echo "process will start at:"
date
echo "++++++++++++++++++++"

cd $SLURM_SUBMIT_DIR

python $baseDir/src/format_qtl_for_coloc.QTLtools.py --qtl_file ${QTLdir}/cis_QTLs/${TISSUE}.cis_all.add_tstat.txt --output $baseDir/input/cis_QTL/${TISSUE}.cis_all.formatted.txt --sample_size $N
' >> ${curr_dir}/submit_formatQTL_${tissue}.slurm &
		wait
		sbatch ${curr_dir}/submit_formatQTL_${tissue}.slurm
	done < $curr_dir/input/sampleN_in_tissues.txt
}


# split qtl by gene
function run_split_qtl_by_gene(){
	dir=/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project/2022-10-10-coloc-smr
	current_dir=`pwd`
	for f in `ls ${dir}/input/cis_QTL/*.cis_all.formatted.txt`
	do
		tissue=`basename $f .cis_all.formatted.txt`
		if [ ! -d "$dir/input/atQTL_byGene/$tissue" ]
		then
			mkdir -p $dir/input/atQTL_byGene/$tissue
		fi

		echo "#!/bin/bash" > ${current_dir}/submit_splitQTL_${tissue}.slurm
		echo "
#SBATCH --job-name=split_$tissue
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=split_${tissue}.err 
#SBATCH --output=split_${tissue}.out " >> ${current_dir}/submit_splitQTL_${tissue}.slurm

		echo $'
echo "process will start at:"
date
echo "+++++++++++++++++++++++++++++"

cd $SLURM_SUBMIT_DIR' >> ${current_dir}/submit_splitQTL_${tissue}.slurm

		echo "
baseDir=$dir
TISSUE=$tissue" >> ${current_dir}/submit_splitQTL_${tissue}.slurm

		echo $'
python $baseDir/src/split_qtl_by_gene.py --qtl_file $baseDir/input/cis_QTL/${TISSUE}.cis_all.formatted.txt --outprefix $baseDir/input/atQTL_byGene/${TISSUE}/

echo "++++++++++++++++++++++++++++++"
echo "process will sleep 10s"
sleep 10
echo "process will end at:"
date
' >> ${current_dir}/submit_splitQTL_${tissue}.slurm &
		wait
		sbatch ${current_dir}/submit_splitQTL_${tissue}.slurm
	done
}



run_atssQuant_by_tissue(){
	tissue=$1
	curr_dir=`pwd`

	echo "#!/bin/bash" > $curr_dir/submit_aTSS_quant.${tissue}.slurm

	echo "
#SBATCH --job-name=$tissue
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=${tissue}.err
#SBATCH --output=${tissue}.out
" >> $curr_dir/submit_aTSS_quant.${tissue}.slurm

	echo $'
echo "process will start at:"
date
echo "++++++++++++++++++++++"

cd $SLURM_SUBMIT_DIR

baseDir=/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project
GTExdir=/lustre/home/xdzou/data/GTEx_BigWig
' >> $curr_dir/submit_aTSS_quant.${tissue}.slurm
	echo "
TISSUE=$tissue" >> $curr_dir/submit_aTSS_quant.${tissue}.slurm

	echo $'
if [ ! -d "$baseDir/output/atss_quant" ]
then
	mkdir -p $baseDir/output/atss_quant
fi

source activate pybigwig_env

python $baseDir/src/Tandem_TSS_calculate.xdzou_v1.py -i $GTExdir/${TISSUE}.bw_list.txt \
	-a $baseDir/input/hg38_transcript_first_exon_anno.txt \
	-t $TISSUE \
	-c 8 \
	-o $baseDir/output/atss_quant &
wait

echo "+++++++++++++++++++++++++"
echo "process will end at:"
date' >> $curr_dir/submit_aTSS_quant.${tissue}.slurm

sbatch $curr_dir/submit_aTSS_quant.${tissue}.slurm
}

#-- run QTLtools permutation and conditional analysis
function run_QTLtool_conditional(){
	curr_dir=`pwd`
	baseDir=/lustre/home/xdzou/2022-08-05-altTSS_QTL-Project
	if [ ! -d "$baseDir/output/QTL_mapping/permutation/fdr" ]
	then
		mkdir -p $baseDir/output/QTL_mapping/permutation/fdr
	fi

	if [ ! -d "$baseDir/output/QTL_mapping/conditional" ]
	then
		mkdir -p $baseDir/output/QTL_mapping/conditional
	fi

	for tissue in `cat ${curr_dir}/input/gtex_tissues_49.txt`
	do
		if [ ! -f "$baseDir/output/QTL_mapping/conditional/${tissue}.conditional.txt" ]
		then
			echo $tissue
			echo "#!/bin/bash" > ${curr_dir}/submit_QTLtools.${tissue}.condition.slurm
			echo "
#SBATCH --job-name=${tissue}_conditional
#SBATCH --partition fat-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=${tissue}_conditional.err
#SBATCH --output=${tissue}_conditional.out

baseDir=$baseDir
VCF=/lustre/home/wychen/2021-09-06-eRNA-QTL/GTEx8/output/vcf_tissue
BIN=/lustre/home/wychen/2021-06-02-Alter-transcript-ends/bin/qtltools/bin
TISSUE=$tissue
" >> ${curr_dir}/submit_QTLtools.${tissue}.condition.slurm
			echo $'
echo "process will start at:"
date
echo "++++++++++++++++++++"

cd $SLURM_SUBMIT_DIR

echo "Start permutation..."
$BIN/QTLtools cis --vcf $VCF/${TISSUE}.maf05.vcf.gz --bed $baseDir/output/QTL_mapping/phenotype/${TISSUE}.phenotype.txt.gz --cov $baseDir/output/QTL_mapping/covariates/${TISSUE}.covariates.txt --permute 1000 --normal --out $baseDir/output/QTL_mapping/permutation/${TISSUE}.permutation.txt &
wait
echo "Permutation end!"

Rscript $baseDir/src/qtltools_runFDR_cis.R $baseDir/output/QTL_mapping/permutation/${TISSUE}.permutation.txt 0.05 $baseDir/output/QTL_mapping/permutation/fdr/${TISSUE} &
wait

echo "Start condition analysis..."
$BIN/QTLtools cis --vcf $VCF/${TISSUE}.maf05.vcf.gz --bed $baseDir/output/QTL_mapping/phenotype/${TISSUE}.phenotype.txt.gz --cov $baseDir/output/QTL_mapping/covariates/${TISSUE}.covariates.txt --mapping  $baseDir/output/QTL_mapping/permutation/fdr/${TISSUE}.thresholds.txt --normal --out $baseDir/output/QTL_mapping/conditional/${TISSUE}.conditional.txt &
wait

echo "+++++++++++++++++++++++++"
echo "process will end at:"
date
' >> ${curr_dir}/submit_QTLtools.${tissue}.condition.slurm &
#			wait
#			sbatch ${curr_dir}/submit_QTLtools.${tissue}.condition.slurm
		fi
 done
}

# -- main
main
