from __future__ import print_function
import sys
import gzip
import pybedtools
import random
from optparse import OptionParser

#parser = OptionParser(usage=usage)


conversion_file = open(sys.argv[1], 'r')

conversion_dict={}
for line in conversion_file:
	#print(line, file=sys.stderr)
	ln = line.strip().split('\t')
	#print(ln[1],file=sys.stderr)
	msp_id = ln[0]
#########
### Include this section for flipping the haplotype assignment
### I assumed that the first haplotype in the vcf file is hap_1, but maybe Ben doesn't read it this way,
### so I can also check overlap with TreeCalls by flipping the haplotype assignment, i.e. sample 6 in TreeCalls is now msp_110_2, not msp_110_1
#	msp_id_diff_list = msp_id.split('_')
#	if msp_id_diff_list[2]=='1':
#		msp_id = msp_id_diff_list[0]+'_'+msp_id_diff_list[1]+'_'+'2'
#	elif msp_id_diff_list[2]=='2':
#		msp_id = msp_id_diff_list[0]+'_'+msp_id_diff_list[1]+'_'+'1'
#########	
	ind_id = ln[1]
	conversion_dict[ind_id]=msp_id

conversion_file.close()
#print(conversion_dict, file=sys.stderr)

Tree_bed_file = gzip.open(sys.argv[2], 'rb')
#Tree_bed_file_mod = open(sys.argv[2]+'.mod' , 'w')

for line in Tree_bed_file:
	ln = line.strip().split('\t')
	chrm = ln[0]
	strt = ln[1]
	end = ln[2]
	ind_id = ln[3]
	msp_id = conversion_dict[ind_id] #msp_1000_1
	msp_id_list = msp_id.split('_')
	msp_id_ch = msp_id_list[0]+'_'+msp_id_list[1]+'_'+chrm+'_'+msp_id_list[2]
	new_ln = msp_id_ch+'\t'+strt+'\t'+end+'\t'+chrm
	print(new_ln,file=sys.stdout)
#	Tree_bed_file_mod.write(new_ln)

Tree_bed_file.close()
#Tree_bed_file_mod.close()
