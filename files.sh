#!/bin/bash
source ~/.bashrc

echo 'move vcfs'
mv *.vcf.* vcfs/

echo 'move window files'
mv *.windowcalc_out.gz RegionFiles/

echo 'make sample_list'
for msp in $(cat Tenn.popfile | grep -v samp | cut -f 1); do echo -e $msp":"1'\n'$msp":"2 >> Tenn.sample_list; done

cat Tenn.sample_list | awk 'BEGIN {OFS="\t"} {print $1, NR-1}' > Tenn.sample_list2 ; mv Tenn.sample_list{2,}

echo 'make chr_list'
ls -C RegionFiles/ | cut -f 3 -d '_' | sort - | uniq - > Tenn.chr_list
