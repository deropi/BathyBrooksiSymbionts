#!/usr/bin/env python2
import sys
import os
import argparse

parser = argparse.ArgumentParser(description = "Neutrality Index and Derived Allele Frequency Spectra estimation", epilog = 'Please, report bugs to dpicazo@ifam.uni-kiel.de')
parser.add_argument("path_to_vcf",metavar="path_to_vcf",type=str,help="Path to the folder containing the vcf files to be analyzed")
parser.add_argument("genome_strains",metavar="genome_strains",type=str,help="Fasta file of the genome analyzed")
parser.add_argument("output",metavar="output",type=str,help="name output folder")
parser.add_argument("strain_clades",metavar="strain_clades",type=str,help="strain clades relational table")
parser.add_argument("ancestor",metavar="ancestor",type=str,help="ancestral strain")
parser.add_argument("good_snps",metavar="good_snps",type=str,help="good snps")
args=parser.parse_args()

################################################# Reading list of non-considered SNPs #######################################
good = {}

with open(args.good_snps, 'r') as remove:
	for line in remove:
		snp = line.strip()
		good[snp] = None


############################################### Genecode
gencode = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}



################################################### Reading and storing the SNPs ###########################################
if not args.path_to_vcf.endswith('/'):
	args.path_to_vcf += '/'
snps_pool = {}	#1st.key:gene, 2nd key:position, value: alt allele from .vcf
frequencies = {}	#1st. key:snp, 2nd key:sample, value:freq from .vcf



states = {}
samples = 0
for file in os.listdir(args.path_to_vcf):
	if file.endswith('.vcf'):
		samples +=1 
		current_sample = file[0]
		with open('%s%s'%(args.path_to_vcf, file), 'r') as vcf:
			for line in vcf:
				if not line.startswith('#'):
					info = line.strip().split('\t')
					gene = info[0]
					position = info[1]
					alt = info[4]
					ref = info[3]
					dp4 = info[7].split(';')[3]
                                        alt_values = int(dp4.split(',')[-1]) + int(dp4.split(',')[-2])
                                        ref_values = int(dp4.split(',')[-4].split('=')[1]) + int(dp4.split(',')[-3])
                                        freq = float(alt_values)/(float(alt_values)+float(ref_values))
					snp_id = gene + '*' + position
					if not snp_id in frequencies:
						frequencies[snp_id] = {}
						frequencies[snp_id][current_sample] = freq
						states[snp_id] = [ref,alt]
					else:
						frequencies[snp_id][current_sample] = freq
					if not gene in snps_pool:
						snps_pool[gene]={}
						snps_pool[gene][position]= alt
					elif not position in snps_pool[gene]:
						snps_pool[gene][position]= alt
					#Here we remove multi-state
					elif snps_pool[gene][position] != alt:
						del snps_pool[gene][position]
################################################# Storing strain sequences ######################################################
genomes = {}

with open(args.genome_strains, 'r') as fasta:
	for line in fasta:
		if line.startswith('>'):
			header= line.strip()[1:]
			strain = header.split('-')[0]
			gene = header.split('-')[1]
			seq = fasta.next().strip()
			if genomes.has_key(strain):
				genomes[strain][gene] = seq
			else:
				genomes[strain] = {}
				genomes[strain][gene] = seq

clades = {}
with open(args.strain_clades, 'r') as clade_file:
	for line in clade_file:
		info = line.strip().split('\t')
		clades[info[0]] = info[1]
################################################## Classifying SNPs into fixed or polymorphic and based on the presence in strains ##############################

snp_categories = {'fixed_syn':0, 'fixed_nonsyn':0, 'polymorphic_invariant_syn':0, 'polymorphic_invariant_nonsyn':0, 'polymorphic_variant_syn':0, 'polymorphic_variant_nonsyn':0}
out = open(args.output, 'w')
out.write('sample\tsnp\tfrequency\tSNP_class\n')

def print_out(frequencies, snp_id, output, type_mutation):
	for sample in frequencies[snp_id]:
		output.write(sample + '\t' + snp_id + '\t' + str(frequencies[snp_id][sample]) + '\t' + type_mutation + '\n')
def print_out_newalt(frequencies, snp_id, output, type_mutation):
        for sample in frequencies[snp_id]:
                output.write(sample + '\t' + snp_id + '\t' + str(1-(frequencies[snp_id][sample])) + '\t' + type_mutation + '\n')
snps_count = 0
not_assigned = set()


## snps_pool structure: 1st-key:gene, 2nd-key:position, value:alt codon
for gene in snps_pool:
	for snp in snps_pool[gene]:
		snp_id = gene + '*' + snp
		status1 = ''
		status2 = ''
		if snp_id in good:
			snps_count += 1
			clades_present = []
			strains_present = []
			snp = int(snp)
			incomplete_groups = 0
			diff_found = 0
			## Check in which stran we have the SNP
			for strain in genomes:
				if int(snp)%3 == 0:
					start = snp-3
					end = snp
					strain_codon = genomes[strain][gene][snp-3:snp]
					snp_codon = genomes[strain][gene][snp-3:snp-1] + snps_pool[gene][str(snp)]
					status1 = genomes[strain][gene][snp-3:snp-1] + states[snp_id][0]
					status2 = genomes[strain][gene][snp-3:snp-1] + states[snp_id][1]
        			elif int(snp)%3==1:
					start = snp -1
					end = snp +2
					strain_codon = genomes[strain][gene][snp-1:snp+2]
					snp_codon = snps_pool[gene][str(snp)] +  genomes[strain][gene][snp:snp+2]
					status1 = states[snp_id][0] + genomes[strain][gene][snp:snp+2]
					status2 = states[snp_id][1] + genomes[strain][gene][snp:snp+2]
        			elif int(snp)%3==2:
					start = snp-2
					end = snp +1
					strain_codon = genomes[strain][gene][snp-2:snp+1]
                			snp_codon = genomes[strain][gene][snp-2] + snps_pool[gene][str(snp)] + genomes[strain][gene][snp]
					status1 = genomes[strain][gene][snp-2] + states[snp_id][0] + genomes[strain][gene][snp]
					status2 = genomes[strain][gene][snp-2] + states[snp_id][1] + genomes[strain][gene][snp]
				if strain_codon == snp_codon:
					clades_present.append(clades[strain])
					strains_present.append(strain)
				else:
					diff_found += 1

################# This part is meant to decide whether an allele is ANCESTRAL or DERIVED
			## If the number of strains with the ALT allele over
			if len(set(clades.keys()) - set(strains_present)) < len(strains_present):
				strains_present = list(set(clades.keys()) - set(strains_present))
				clades_present = []
				for element in strains_present:
					clades_present.append(clades[element])
				for element in set(clades_present):
					if clades_present.count(element) < clades.values().count(element):
						incomplete_groups +=1
				if len(set(clades_present)) ==1 and  incomplete_groups == 0:
					print_out_newalt(frequencies, snp_id, out, clades_present[0])
				elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
					print_out_newalt(frequencies, snp_id, out, 'multiple_groups')
			elif len(set(clades.keys()) - set(strains_present)) > len(strains_present): #len(strains_present) > len(set(clades.keys()) - set(strains_present)):
                                for element in set(clades_present):
                                        if clades_present.count(element) < clades.values().count(element):
                                                incomplete_groups +=1
                                if len(set(clades_present)) ==1 and  incomplete_groups == 0:
                                        print_out(frequencies, snp_id, out, clades_present[0])
                                elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
                                        print_out(frequencies, snp_id, out, 'multiple_groups')
			elif len(strains_present) > 0:
				if args.ancestor in strains_present:
					strains_present = list(set(clades.keys()) - set(strains_present))
                                	clades_present = []
                                	for element in strains_present:
                                        	clades_present.append(clades[element])
					for element in set(clades_present):
	                                        if clades_present.count(element) < clades.values().count(element):
        	                                        incomplete_groups +=1
					if len(set(clades_present)) ==1 and  incomplete_groups == 0:
                                        	print_out_newalt(frequencies, snp_id, out, clades_present[0])
                                	elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
                                        	print_out_newalt(frequencies, snp_id, out, 'multiple_groups')
				else:
					for element in set(clades_present):
	                                        if clades_present.count(element) < clades.values().count(element):
        	                                        incomplete_groups +=1
					if len(set(clades_present)) ==1 and  incomplete_groups == 0:
	                                        print_out(frequencies, snp_id, out, clades_present[0])
        	                        elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
                	                        print_out(frequencies, snp_id, out, 'multiple_groups')
			if len(strains_present) == 0:
                                snp_type = 'polymorphic_invariant'
                                print_out(frequencies, snp_id, out, 'polymorphic_invariant')
			elif incomplete_groups == 0 and len(strains_present) > 0:
				snp_type = 'fixed'
			elif incomplete_groups > 0 and len(strains_present) > 0:
				snp_type = 'polymorphic_variant'
                                print_out(frequencies, snp_id, out, 'polymorphic_variant')
			if gencode[status1] != gencode[status2]:
				snp_type += '_nonsyn'
			else:
				snp_type += '_syn'
              		snp_categories[snp_type] +=1
out.close()
print snps_count
print snp_categories
