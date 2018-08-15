#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=20G
#SBATCH --qos=1day
#SBATCH --time=23:00:00
#SBATCH --output=/scratch/tmp/abwolf/Sstar/null/slrun.Sstar.Tenn.null.sh.%A_%a.o
source ~/.bashrc

date
echo $SLURM_JOB_NAME
echo $SLURM_ARRAY_TASK_ID
echo ''

mdl=$( echo Tenn)
seed=$( echo ${SLURM_ARRAY_TASK_ID}${RANDOM} )
n1=$( echo 0.0)
n2=$( echo 0.0)
eur=$( echo 1006 )   #1006
asn=$( echo 1008 )   #1008 ; 2040
len=$( awk 'BEGIN {print 10e6}' )

tag=$(echo "$mdl"_nonAfr_"$seed"_n1_"$n1"_n2_"$n2")
dir=~/SimulatedDemographic/msprime/
sstardir=~/SimulatedDemographic/Sstar/

echo mdl: $mdl	seed: $seed n1: $n1 n2: $n2 eur: $eur asn: $asn		tag: $tag
echo **RUN MSPRIME SIMULATION AND OUTPUT VCF**
python $dir/bin/msprime.Admixture_Simulate.py \
			-p nonAfr \
			-o $mdl \
			-s $seed \
			-i 2 \
			-n $n1 \
			-d $n2 \
			-e $eur \
			-a $asn \
			-r 216 \
			-c vcf \
			-l $len

echo fin simulation
date


echo **MODIFY VCF , BGZIP , TABIX**
echo change chr 
zcat $tag.vcf.gz | awk 'BEGIN {OFS="\t"} /^#/{print$0} !/^#/{$1="'$seed'" ; print $0}' | bgzip -c > $tag.mod.vcf.gz

echo tabix
tabix -p vcf $tag.mod.vcf.gz

echo remove original
rm $tag.vcf.gz
echo fin
date


#echo **REMOVE SINGLETON-DOUBLETON-TRIPLETON SITES FROM VCF**
#vcftools --gzvcf $tag.mod.vcf.gz \
#	--mac 4 \
#	--recode --stdout \
#	| bgzip -c \
#	> $tag.mod.filtered.vcf.gz
#
#echo rename filtered file
#mv $tag.mod.filtered.vcf.gz $tag.mod.vcf.gz
#tabix -p vcf $tag.mod.vcf.gz



echo **EXTRACT ARCHAIC DATA FROM VCF**
vcftools --gzvcf $tag.mod.vcf.gz \
	--keep $sstardir/bin/vcf_keep_archaic.txt \
	--recode --stdout \
	| bgzip -c \
	> $tag.mod.ARCHAIC.vcf.gz
tabix -p vcf $tag.mod.ARCHAIC.vcf.gz

echo fin
date

echo **CALCULATE WINDOWED SSTAR SCORES AND MATCH PCT -- updated Sstar version**
~/software/anaconda2/bin/python $sstardir/s_star_git/bin/windowed_calculations.py \
	--vcf-has-illumina-chrnums \
	-vcfz $tag.mod.vcf.gz \
	-indf ./$mdl.popfile \
	-target-pops EUR ASN \
	-ref-pops AFR \
	--archaic-vcf $tag.mod.ARCHAIC.vcf.gz \
	-p 10 \
	-s-star \
	-winlen 50000 \
	-winstep 10000 \
	-no-pvalues \
	-range 0 $len \
	| gzip -c - \
	> $tag.windowcalc_out.gz



echo ''
echo fin Sstar calc
date


#echo **CONVERT REGIONS FILE TO INDFILES**
##
##mkdir ./RegionFiles/ ./IndFiles/ ./SstarSigFiles/
##
##zcat $tag.windowcalc_out.gz | awk 'BEGIN {OFS="\t"} {if(NR!=1) print $7}' | sort - | uniq - > $mdl.sample_list
##
#gzip -dc $tag.windowcalc_out.gz > RegionFiles/$tag.windowcalc_out
#
#
#perl $sstardir/bin/windowcalc_single_file_to_indfiles.pl \
#	. \
#	./$mdl.sample_list \
#	$tag \
#
#echo ''
#echo fin
#date
#
#echo **CALCULATE SSTAR P-VALUES FROM INDFILES**
#
#for asamp in $(cat ./$mdl.sample_list); do
#	echo $asamp
#	Rscript $sstardir/bin/compute_pvalues_from_indfiles.r \
#		$tag \
#		$sstardir/bin/Tenn.neand_callable_bases_chr1to100_10Mb \
#		$sstardir/bin/Tenn.recomb_rate_per_50kb_window_chr1to100_10Mb \
#		. \
#		glm \
#		$sstardir/bin/gravel_asn_scale_60k.model.Rdata \
#		$asamp
#		
#		# tag for file
#		# bed file of Neand callable bases
#			##chr winstart  winend  callable SNPs_per_50kb_window  <-- This column does not matter, it gets dropped when table is loaded into script
#			##10	0	50000	50000	650
#			##10	10000	60000	50000	650
#			##10	20000	70000	50000	650
#			##10	30000	80000	50000	650
#		#recomb rate file
#			##chr	winstart recomb_rate_cM_per_MB
#			##10	0	1
#			##10	50000	1
#			##10	100000	1
#		#output directory
#		#model type, either glm or ecdf
#		#sstar_null_model
#		#sample_name
#done
#echo fin
#date
#
#echo **CALCULATE Q-VALUE THRESHOLDS FOR SSTAR AND MATCH PCT**
#
#Rscript $sstardir/bin/calculate_qvalue_thresholds.r \
#	$tag \
#	./$mdl.sample_list \
#	. \
#	5 \
#	
#	# tag for file
#	# sample list --> only those that have a Sstarsig file with data
#	# output directory
#	# minimum n_Sstar snps in window
#
echo FIN
date
