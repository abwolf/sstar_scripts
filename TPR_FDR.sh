#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=10G
#SBATCH --qos=1day
#SBATCH --time=22:00:00
source ~/.bashrc

date

admix=$( echo n1_0.1_mAfB_0.0_mBAf_0.0_mAfEu_0.0_mEuAf_0.0 )

for i in $( cat Tenn.chr_list ); do
	echo $i
	sstarf=$( echo bedfiles/Tenn_nonAfr_"$i"_"$admix".sstar_sig_0.01.match_sig_N_0.05.isc_0.bed )
	TreeCallf=$( echo TreeCalls/Tenn_nonAfr_"$i"_"$admix".bed.mod.merged )

	gzip -dv $TreeCallf.gz

	cat $sstarf | grep -v msp_ID | sort-bed - \
	| bedmap --ec --delim '\t' --echo --bases-uniq-f - $TreeCallf \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1) FD+=1 ; if($4>=0.1) TD+=1} END {print "chr: "'$i', "TD: "TD, "FD: "FD, "FDR: "FD/(TD+FD)}' \
	>> TPR_FDR.txt

    cat $sstarf | grep -v msp_ID | sort-bed - \
    | bedmap --ec --delim '\t' --echo --bases-uniq-f $TreeCallf - \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1) FN+=1 ; if($4>=0.1) TP+=1} END {print "chr: "'$i', "TP: "TP, "FN: "FN, "TPR: "TP/(TP+FN)}' \
	>> TPR_FDR.txt


	cat $sstarf | grep -v msp_ID | sort-bed - \
	| bedmap --ec --delim '\t' --echo --bases-uniq-f - $TreeCallf \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1 && $3-$2>=30000) FD+=1 ; if($4>=0.1 && $3-$2>=30000) TD+=1} END {print "chr: "'$i', "TD: "TD, "FD: "FD, "FDR: "FD/(TD+FD)}' \
	>> TPR_FDR.30kb.txt

    cat $sstarf | grep -v msp_ID | sort-bed - \
    | bedmap --ec --delim '\t' --echo --bases-uniq-f $TreeCallf - \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1 && $3-$2>=30000) FN+=1 ; if($4>=0.1 && $3-$2>=30000) TP+=1} END {print "chr: "'$i', "TP: "TP, "FN: "FN, "TPR: "TP/(TP+FN)}' \
	>> TPR_FDR.30kb.txt

	gzip -v $TreeCallf
done

date
echo FIN
