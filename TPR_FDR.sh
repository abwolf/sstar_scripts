#!/bin/bash
#SBATCH --get-user-env
#SBATCH --mem=10G
#SBATCH --qos=1day
#SBATCH --time=22:00:00
source ~/.bashrc

date

admix=$admix

#admix=$1

for i in $( cat $admix/Tenn.chr_list ); do
	echo $i
	sstarf=$( echo $admix/bedfiles/Tenn_nonAfr_"$i"_"$admix".sstar_sig_0.01.match_sig_N_0.05.isc_0.bed )
	TreeCallf=$( echo $admix/TreeCalls/Tenn_nonAfr_"$i"_"$admix".bed.mod )

	zcat $TreeCallf.gz | sort-bed - > $TreeCallf

	cat $TreeCallf | awk 'BEGIN {OFS="\t"} {if($3-$2>30000) print $0}' > $TreeCallf.30kb

	cat $sstarf | grep -v msp_ID | sort-bed - \
	| bedmap --ec --delim '\t' --echo --bases-uniq-f - $TreeCallf \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1) FD+=1 ; if($4>=0.1) TD+=1} END {print "chr: "'$i', "TD: "TD, "FD: "FD, "FDR: "FD/(TD+FD)}' \
	>> $admix/TPR_FDR.txt

    cat $sstarf | grep -v msp_ID | sort-bed - \
    | bedmap --ec --delim '\t' --echo --bases-uniq-f $TreeCallf - \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1) FN+=1 ; if($4>=0.1) TP+=1} END {print "chr: "'$i', "TP: "TP, "FN: "FN, "TPR: "TP/(TP+FN)}' \
	>> $admix/TPR_FDR.txt


	cat $sstarf | grep -v msp_ID | sort-bed - \
	| bedmap --ec --delim '\t' --echo --bases-uniq-f - $TreeCallf.30kb \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1) FD+=1 ; if($4>=0.1) TD+=1} END {print "chr: "'$i', "TD: "TD, "FD: "FD, "FDR: "FD/(TD+FD)}' \
	>> $admix/TPR_FDR.30kb.txt

    cat $sstarf | grep -v msp_ID | sort-bed - \
    | bedmap --ec --delim '\t' --echo --bases-uniq-f $TreeCallf.30kb - \
	| awk 'BEGIN {OFS="\t"} {if($4<0.1) FN+=1 ; if($4>=0.1) TP+=1} END {print "chr: "'$i', "TP: "TP, "FN: "FN, "TPR: "TP/(TP+FN)}' \
	>> $admix/TPR_FDR.30kb.txt

	rm $TreeCallf $TreeCallf.30kb
done

date
echo FIN
