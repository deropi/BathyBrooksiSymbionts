#!/usr/bin/env python2
import sys
import os

path_to_vcf = sys.argv[1]
genome_strains = sys.argv[2]
good_snps = sys.argv[3]
symbiont = sys.argv[4]
output = sys.argv[5]


################################################ New sample naming ########################################################
annotations = {}

with open('/home/dpicazo/Population-structure/new_annot.txt', 'r') as table:
        for line in table:
                if not line.startswith('old'):
                        info = line.strip().split('\t')
                        annotations[info[0]] = info[1]

################################################# Reading list of non-considered SNPs #######################################
good = {}

with open(good_snps, 'r') as remove:
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
if not path_to_vcf.endswith('/'):
	path_to_vcf += '/'
snps_pool = {}	#1st.key:gene, 2nd key:position, value: alt allele from .vcf
frequencies = {}	#1st. key:snp, 2nd key:sample, value:freq from .vcf



states = {}
samples = 0
for file in os.listdir(path_to_vcf):
	if file.endswith('.vcf'):
		samples +=1 
		current_sample = file[0]
		with open('%s%s'%(path_to_vcf, file), 'r') as vcf:
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

with open(genome_strains, 'r') as fasta:
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
if symbiont == 'MOX':
	clades = {'strain_2_aligned': 'group1', 'strain_1_aligned': 'group2', 'strain_4_aligned': 'group1','strain_3_aligned': 'group2', 'strain_5_aligned': 'group2', 'strain_6_aligned': 'group1'}
	ancestor = 'strain_1_aligned'
else:
	clades = {'strain_2_aligned': 'group1', 'strain_1_aligned': 'group1', 'strain_4_aligned': 'group1', 'strain_10_aligned': 'group1', 'strain_3_aligned': 'group3', 'strain_7_aligned': 'group3', 'strain_9_aligned': 'group3',  'strain_5_aligned': 'group2', 'strain_6_aligned': 'group2', 'strain_11_aligned': 'group2', 'strain_8_aligned': 'group4'}
	ancestor = 'strain_8_aligned'
################################################## Classifying SNPs into fixed or polymorphic and based on the presence in strains ##############################

snp_categories = {'fixed_syn':0, 'fixed_nonsyn':0, 'polymorphic_invariant_syn':0, 'polymorphic_invariant_nonsyn':0, 'polymorphic_variant_syn':0, 'polymorphic_variant_nonsyn':0}
out = open(output, 'w')
out.write('sample\tsnp\tfrequency\tSNP_class\n')

def print_out(frequencies, snp_id, output, type_mutation, annotations):
	for sample in frequencies[snp_id]:
		output.write(annotations[sample] + '\t' + snp_id + '\t' + str(frequencies[snp_id][sample]) + '\t' + type_mutation + '\n')
def print_out_newalt(frequencies, snp_id, output, type_mutation, annotations):
        for sample in frequencies[snp_id]:
                output.write(annotations[sample] + '\t' + snp_id + '\t' + str(1-(frequencies[snp_id][sample])) + '\t' + type_mutation + '\n')
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
			## This variable stores the number of clades that are incomplete for a specific SNP
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

			### WHAT?? Should not even be working...there are more samples than strains...stupid Devani...
			if diff_found == samples:#diff_found == 0 or diff_found == samples:
				for sample in frequencies[snp_id]:
					frequencies[snp_id][sample] = 1-freq

################# This part is meant to decide whether an allele is ANCESTRAL or DERIVED
			## If the number of strains with the ALT allele exceeds the number of strains with the reference, then the ALT becomes the reference
			elif len(set(clades.keys()) - set(strains_present)) < len(strains_present):
				strains_present = list(set(clades.keys()) - set(strains_present))
				clades_present = []
				for element in strains_present:
					clades_present.append(clades[element])
				for element in set(clades_present):
					if clades_present.count(element) < clades.values().count(element):
						incomplete_groups +=1
				## If there is no incomplete clade, we output the info about the snp for each sample
				if len(set(clades_present)) ==1 and  incomplete_groups == 0:
					print_out_newalt(frequencies, snp_id, out, clades_present[0], annotations)
				elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
					print_out_newalt(frequencies, snp_id, out, 'multiple_groups', annotations)
			## If the number of strains without the SNP exceeds the number of strains with the SNP, then we proceed printing without re-defining the ALT
			elif len(set(clades.keys()) - set(strains_present)) > len(strains_present): #len(strains_present) > len(set(clades.keys()) - set(strains_present)):
                                for element in set(clades_present):
                                        if clades_present.count(element) < clades.values().count(element):
                                                incomplete_groups +=1
                                if len(set(clades_present)) ==1 and  incomplete_groups == 0:
                                        print_out(frequencies, snp_id, out, clades_present[0], annotations)
                                elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
                                        print_out(frequencies, snp_id, out, 'multiple_groups', annotations)
			## If the number of strains with the SNPs equals the number of strains without: 
			elif len(strains_present) > 0:
				## If the defined ancestor has the SNP, then the SNP turns DERIVED, and we output info the same way
				if ancestor in strains_present:
					strains_present = list(set(clades.keys()) - set(strains_present))
                                	clades_present = []
                                	for element in strains_present:
                                        	clades_present.append(clades[element])
					for element in set(clades_present):
	                                        if clades_present.count(element) < clades.values().count(element):
        	                                        incomplete_groups +=1
					if len(set(clades_present)) ==1 and  incomplete_groups == 0:
                                        	print_out_newalt(frequencies, snp_id, out, clades_present[0], annotations)
                                	elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
                                        	print_out_newalt(frequencies, snp_id, out, 'multiple_groups', annotations)
				## If not, then output without re-defining the ancestor
				else:
					for element in set(clades_present):
	                                        if clades_present.count(element) < clades.values().count(element):
        	                                        incomplete_groups +=1
					if len(set(clades_present)) ==1 and  incomplete_groups == 0:
	                                        print_out(frequencies, snp_id, out, clades_present[0], annotations)
        	                        elif len(set(clades_present)) > 1 and  incomplete_groups == 0:
                	                        print_out(frequencies, snp_id, out, 'multiple_groups', annotations)
			## If the SNP is not present in any strain, then it is defined as INVARIANT
#################################################### In adition, here we defined the SNP type...
			if len(strains_present) == 0:
                                snp_type = 'polymorphic_invariant'
                                print_out(frequencies, snp_id, out, 'polymorphic_invariant', annotations)
			## If it is not invariant and groups are complete for the presence of the SNP, this will be a fixed SNP
			elif incomplete_groups == 0 and len(strains_present) > 0:
				snp_type = 'fixed'
				print strains_present
			## If it is partially present among clusters, then we have a polymorphic variant
			elif incomplete_groups > 0 and len(strains_present) > 0:
				snp_type = 'polymorphic_variant'
                                print_out(frequencies, snp_id, out, 'polymorphic_variant', annotations)

			######## HERE, we want in addition include to the snp_type category, whether the snp is synonymous or non_synonymous
			if gencode[status1] != gencode[status2]:
				print gencode[status1], gencode[status2]
				snp_type += '_nonsyn'
			else:
				print gencode[status1], gencode[status2]
				snp_type += '_syn'
              		snp_categories[snp_type] +=1
out.close()
print snps_count
print snp_categories

