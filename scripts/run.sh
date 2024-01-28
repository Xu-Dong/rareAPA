#!/bin/bash

# --------------------- Functions ----------------------
# main function

main(){
	#>> apa outliers pipe
#	run_inds_extractor_apa
#	run_inds_extractor_covariate
#	run_scaling
#	run_gather_matrix
#	run_call_outliers gtex_aOutlier_v8_normalized_pdui.peer..txt apa
#	run_call_singlez_outliers apa

}



# keep only EA indis in PDUI table
function run_inds_extractor_apa(){
	indir=${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/input/peer_pdui
	current_dir=`pwd`
	for f in `ls ${indir}/*.pdui.peer.txt`
	do
		fname=`basename $f`
		prefix=${fname%.pdui.peer.txt}

		echo "
#SBATCH --job-name=INDS_Extraction_$prefix
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=INDS_Extraction_${prefix}.err
#SBATCH --error=INDS_Extraction_${prefix}.out


##################################
indir=$indir
fname=$fname" > ${current_dir}/submit_onlyEA_${prefix}.sh

		echo $'
cd $SLURM_SUBMIT_DIR

Rscript ${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/src/extract_european_individuals.apa.R ${indir}/$fname' >> ${current_dir}/submit_onlyEA_${prefix}.sh &
		wait
		sbatch ${current_dir}/submit_onlyEA_${prefix}.sh &
		wait
	done
}


# keep only EA indis in covariates table
function run_inds_extractor_covariate(){
	indir=${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/input/GTExV8_Covariates
	current_dir=`pwd`
	for f in `ls ${indir}/*.v8.covariates.txt`
	do
		fname=`basename $f`
		prefix=${fname%.v8.covariates.txt}
		echo "
#SBATCH --job-name=INDS_Extraction_$prefix
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=INDS_Extraction_${prefix}.err
#SBATCH --error=INDS_Extraction_${prefix}.out

##################################
indir=$indir
fname=$fname" > ${current_dir}/submit_onlyEA_${prefix}.sh

		echo $'
Rscript ${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/src/extract_european_inds.covariates.R ${indir}/$fname' >> ${current_dir}/submit_onlyEA_${prefix}.sh &
		wait
		sbatch ${current_dir}/submit_onlyEA_${prefix}.sh &
		wait
	done
}

# get Z-transformed PDUI, invoke: scaling_PDUI.R Tissue_name
function run_scaling(){
	basedir=${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/
	current_dir=`pwd`
	for tissue in `cat ${basedir}/input/gtex_v8.eQTL_tissues_49.txt`
	do
		echo "#!/bin/bash" > ${current_dir}/submit_scaling_${tissue}.sh
		echo "
#SBATCH --job-name=scaling_$tissue
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=scaling_${tissue}.err
#SBATCH --output=scaling_${tissue}.out



################################
basedir=$basedir
TISSUE=$tissue" >> ${current_dir}/submit_scaling_${tissue}.sh
		echo $'

cd $SLURM_SUBMIT_DIR

echo "process start at:"
date

Rscript ${basedir}/src/scaling_PDUI.all_samples.R $TISSUE
echo "process end at:"
date
' >>${current_dir}/submit_scaling_${tissue}.sh 
		sbatch ${current_dir}/submit_scaling_${tissue}.sh
	done
}

# gathering PDUI matrix of all tissues together
function run_gather_matrix(){
	basedir=${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA
	/lustre/home/xdzou/src/anaconda2/bin/python ${basedir}/src/gather_filter_normalized_pdui_peerless.py \
	${basedir}/input/gtex_v8.eQTL_tissues_49.txt \
	${basedir}/input/individuals.only_EA.txt \
	.pdui.peer.ztrans.txt \
	${basedir}/output/gtex_aOutlier_v8_normalized_pdui.peer..txt
}

# call multi-tissue aOutliers
function run_call_outliers(){
	input_file=$1
	molecular_type=$2
	basedir=${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA
	curr_dir=`pwd`
	echo "#!/bin/bash" > ${curr_dir}/submit_call_outlier_${molecular_type}.sh
	echo "
#SBATCH --job-name=call_outliers_${molecular_type}
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --error=call_outliers_${molecular_type}.err
#SBATCH --output=call_outliers_${molecular_type}.out

######################" >> ${curr_dir}/submit_call_outlier_${molecular_type}.sh

	echo "
INPUT=$input_file
MOL=$molecular_type
BASEDIR=$basedir" >> ${curr_dir}/submit_call_outlier_${molecular_type}.sh

	echo $'
cd $SLURM_SUBMIT_DIR

Rscript ${BASEDIR}/src/call_outliers_medz_peerless.xdzou.v3.R $INPUT $MOL' >> ${curr_dir}/submit_call_outlier_${molecular_type}.sh 
	sbatch ${curr_dir}/submit_call_outlier_${molecular_type}.sh 
}

# call single-tissue aOutliers
function run_call_singlez_outliers(){
	molecular_type=$1
	bindir=${HOME}/2021-01-01-RareVar-aQTL-Project/2021-07-28-call-outliers-onlyEA/src
	curr_dir=`pwd`
	echo "#PBS -N call_outliers_singlez_${molecular_type}
#PBS -q fat-1
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00
#PBS -j oe

#PBS -V 
#PBS -S /bin/bash
###################
BIN=${bindir}" > ${curr_dir}/submit_call_singlez_${molecular_type}.sh
	if [ $molecular_type == "expression" ]
	then
		echo $'
python ${BIN}/call_outliers_single_tissue.py' >> ${curr_dir}/submit_call_singlez_${molecular_type}.sh &
	else
		echo $'
python ${BIN}/call_outliers_single_tissue.apa.py' >> ${curr_dir}/submit_call_singlez_${molecular_type}.sh &
	fi

	wait
	qsub ${curr_dir}/submit_call_singlez_${molecular_type}.sh &
	wait
}


# ------------- Main ---------------------
main
