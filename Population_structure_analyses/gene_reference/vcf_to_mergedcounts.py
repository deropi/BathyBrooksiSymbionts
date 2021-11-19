#!/usr/bin/env python2
import sys
import os
import re
import argparse

######### Argparse
parser = argparse.ArgumentParser(description = "This script merges per sample SNV data from .vcf files into counts table required for further population structure analyses. The standard format is Lofreq output containing DP4 information", epilog = 'Please, report bugs to dpicazo@ifam.uni-kiel.de')
parser.add_argument("path",metavar="path",type=str,help="Folder containing .vcf files to be parsed. Also, output will be written in the same directory")
args=parser.parse_args()
##########################
if not args.path.endswith('/'):
	args.path += '/'
all_files = os.listdir(args.path)
files = [i for i in all_files if i.endswith('.vcf')]
#########################
all_positions = {}
counts = {}

def read_file(vcf, all_positions, path, counts):
        with open('%s%s' %(path, vcf), 'r') as input_vcf:
		counts[vcf] = {}
#{'A':0, 'C':0, 'G':0, 'T':0}
		for line in input_vcf:
			if not line.startswith('#') and not re.search('INDEL', line):
				info = line.strip().split('\t')
				#print line
				chr_pos = info[0] + '-' + info[1]
				counts[vcf][chr_pos] = {'A':0, 'C':0, 'G':0, 'T':0}
				ref_allele = info[3]
				alt_allele = info[4]
				dp4 = info[7].split(';')[3]
				alt_values = int(dp4.split(',')[-1]) + int(dp4.split(',')[-2])
				ref_values = int(dp4.split(',')[-4].split('=')[1]) + int(dp4.split(',')[-3])
				all_positions[chr_pos] = info[3]
				counts[vcf][chr_pos][alt_allele] += alt_values
				counts[vcf][chr_pos][ref_allele] += ref_values
				
###########################
if len(files) > 0:	
	for vcf in files:
		read_file(vcf, all_positions, args.path, counts)

else:
	print 'Oops!...no .vcf file found in the directory'


#print counts

nt_order = 'ACGT'
for sample in counts:
	sample_name = sample.split('.')[0]
	with open('%s%s_SNPs.list' %(args.path,sample_name), 'w') as out:
		out.write('CHROM\tPos\tRef\tA\tC\tG\tT\n')
		for snp in all_positions.keys():
			snp_info = snp.split('-')
			if counts[sample].has_key(snp):
				out.write(snp_info[0] + '\t'+ snp_info[1] + '\t' + all_positions[snp])
				for nt in nt_order:
					out.write('\t' + str(counts[sample][snp][nt]))
				out.write('\n')
			else:
				out.write(snp_info[0] + '\t'+ snp_info[1] + '\t' + all_positions[snp])
				for nt in nt_order:
					if nt == all_positions[snp]:
						out.write('\t1')
					else:
						out.write('\t0')
				out.write('\n')
	
