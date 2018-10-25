#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=50G
#SBATCH --qos=1day
#SBATCH --time=23:00:00
#SBATCH --output=./slrun.SstarECDFpvalueCalculation.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

sstardir=~/SimulatedDemographic/Sstar/
dir=/Genomics/akeylab/abwolf/SimulatedDemographic/Sstar/deserts
mdl=$( echo Tenn )
admix=$( echo n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0001 )

cmd=$( echo " Rscript $sstardir/bin/SstarECDFpvalueCalculation.R \n
	--inputdir $dir/$mdl/ \n
	--outputdir $dir/$mdl/$admix/bedfiles/  \n
	--mdl '$mdl'_nonAfr \n
	--admix_dir /$admix/ \n
	--admix_tag $admix \n
	--max_chrm_admix 1000 \n
	--ecdf $dir/$mdl/null/SstarECDF_maxchrm_9500.RData.gz \n
	--filter " )

# cmd=$( echo " Rscript $sstardir/bin/SstarECDFpvalueCalculation.R \n
# 	--inputdir $dir/$mdl/ \n
# 	--null_dir /null/ \n
# 	--mdl '$mdl'_nonAfr \n
# 	--null_tag n1_0.0_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0 \n
# 	--max_chrm_null 9500 \n
# 	--ecdf_only " )

echo -e $cmd

eval $( echo -e $cmd )

echo FIN
date
