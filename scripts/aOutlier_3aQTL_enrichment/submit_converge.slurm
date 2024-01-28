#!/bin/bash 

#SBATCH --job-name=converge
#SBATCH --partition cu-1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=converge.err
#SBATCH --output=converge.out 

echo "process will start at:"
date
echo "+++++++++++++++++++++++++++++"

cd $SLURM_SUBMIT_DIR

baseDir=/lustre/home/xdzou/2021-01-01-RareVar-aQTL-Project/2023-01-15-aQTL_aOutlier
for TISSUE in `cat $baseDir/input/gtex_49_tissues.txt`:
do
	echo $TISSUE
	Rscript $baseDir/src/check_aOutliers_aQTL_bootstrap.R -t $TISSUE &
	wait
	Rscript $baseDir/src/check_aOutliers_aQTL_bootstrap.background.R -t $TISSUE &
	wait
done

echo "++++++++++++++++++++++++++++++"
echo "process end at:"
date

