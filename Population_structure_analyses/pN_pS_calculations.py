#!/usr/bin/env python2
import sys
import argparse
import os

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

parser = argparse.ArgumentParser(description = "sample-specific pN/pS estimation from vcf files and genome fasta file", epilog = 'Please, report bugs to dpicazo@ifam.uni-kiel.de')
parser.add_argument("path_to_vcf",metavar="path_to_vcf",type=str,help="Path to the folder containing the vcf files to be analyzed")
parser.add_argument("genome",metavar="genome",type=str,help="Fasta file of the genome analyzed")
parser.add_argument("out_name",metavar="out_name",type=str,help="name output folder")
#parser.add_argument("freq",metavar="freq",type=str,help="frequency threshold for fixed snps")
parser.add_argument("good",metavar="good",type=str,help="good snps")
args=parser.parse_args()

# We store a list of SNPs that we want to analyze
good_ones = []
with open(args.good, 'r') as goods:
	for line in goods:
		snp = line.strip()
		good_ones.append(snp)


########################## THIS PIECE OF CODE IS MEANT TO PREPARE THE DATA STRUCTURE

# Here we create 'possible_mutations'dict, which is saving ALL possible syn and non-syn changes for each codon in the genecode
possible_mutations = {}

nts = 'AGCT'

for codon in gencode:
        codon_list = []
        possible_mutations[codon] = {}
        possible_mutations[codon]['syn']= 0
        possible_mutations[codon]['non_syn']= 0
        codon_list = [e[0] for e in codon]
        for i in range(0, len(codon_list)):
                for nt in nts:
                        new_codon_list =[e[0] for e in codon_list]
                        #Here we have mutated our codon
			new_codon_list[i] = nt
			#Here we make sure that we only count mutations
                        if new_codon_list != codon_list:
                                new_codon = ''.join(new_codon_list)
                                if gencode[codon] == gencode[new_codon]:
                                        possible_mutations[codon]['syn'] +=1
                                else:
                                        possible_mutations[codon]['non_syn']+= 1

#This dictionary is meant to store the two states that we can find for each of the mutated positions
# 1st_key: 'codon_position*position in the codon';  2nd.key: state; value: codon
genome_positions = {}
#Here, we store the REFERENCE GENOME
genome = {}
#Here, we store the possible mutations for each position in our genome, for genes or genomewise
genome_wise_possible = {'syn':0, 'non_syn':0}
per_gene_possible = {}
with open(args.genome, 'r') as fasta:
	for line in fasta:
		if line.startswith('>'):
			header = line[1:].strip()
			seq = fasta.next().strip()
			genome[header] = seq
			per_gene_possible[header] = {}
			per_gene_possible[header]['syn'] = 0
			per_gene_possible[header]['non_syn'] =0
			position = 1
			genome_positions[header] = {}
			for i in range(0, len(seq), 3):
				codon = seq[i:i+3]
				genome_wise_possible['syn'] += possible_mutations[codon]['syn']
				genome_wise_possible['non_syn'] += possible_mutations[codon]['non_syn']
				per_gene_possible[header]['syn'] += possible_mutations[codon]['syn']
				per_gene_possible[header]['non_syn'] += possible_mutations[codon]['non_syn']
				for x in range(0, 3):
					genome_positions[header][str(position) + '*' + str(x+1)] = {}
					genome_positions[header][str(position) + '*' + str(x+1)]['state1'] = codon
				position +=1



########################################################### FROM NOW ON WE ARE READING OUR DATA


multistate_alleles= {}
sample_data = {}


samples = []



#Ths function is retreeving the codon from the snp information
def snp_to_codon(locus, gene, alt_snp, genome_positions):
	codon_pos = ((int(locus)-1)/3)+1
	if int(locus)%3 == 0:
		reference_codon = genome_positions[gene][str(codon_pos) + '*3']['state1']
	else:
		reference_codon = genome_positions[gene][str(codon_pos) + '*' + str(int(locus)%3)]['state1']
	reference_codon_list = list(reference_codon)
	if int(locus)%3==1:
		reference_codon_list[0] = alt_snp
	elif int(locus)%3==2:
		reference_codon_list[1] = alt_snp
	else:
		reference_codon_list[2] = alt_snp
	return ''.join(reference_codon_list)


# This one is obvious...
def read_vcf(path_to_vcf, genome_positions, multistate_alleles, file_name, sample_data, samples, good_ones):
	with open('%s%s'%(path_to_vcf,file_name), 'r') as vcf:
		sample = file_name[0]
		samples.append(sample)
		for line in vcf:
			if not line.startswith('#'):
				info = line.strip().split('\t')
				#Here we get the codon number
				codon_position = ((int(info[1])-1)/3)+1
				tag = info[0] + '*' + info[1]
				#The conditions here are for codon positions
				if tag in good_ones:
					if int(info[1])%3 ==0:
						snp = str(codon_position) + '*3'
					else:
						snp = str(codon_position) + '*' + str(int(info[1])%3)
					alt = info[4]
					#Here we are mutating the codon
					alt_codon = snp_to_codon(info[1], info[0], alt, genome_positions)
					## These two next conditions are just check-points
					if not genome_positions[info[0]].has_key(snp):
						print 'weird...'
					if genome_positions[info[0]][snp].has_key('state2') and genome_positions[info[0]][snp]['state2'] != alt_codon and genome_positions[info[0]][snp]['state1'] != alt_codon:
						multistate_alleles[str(info[0]) + '*' + snp] = None
					elif not genome_positions[info[0]][snp].has_key('state2'):
						genome_positions[info[0]][snp]['state2'] = alt_codon
					if not sample_data.has_key(info[0] + '*' + snp):
						sample_data[info[0] + '*' +snp] = {}
						sample_data[info[0] + '*' +snp][sample] = {}
					else:
						sample_data[info[0] + '*' +snp][sample] = {}
					dp4 = info[7].split(';')[3]
					alt_values = int(dp4.split(',')[-1]) + int(dp4.split(',')[-2])
					ref_values = int(dp4.split(',')[-4].split('=')[1]) + int(dp4.split(',')[-3])
					if ref_values > 0:
						frequency_alt =  float(alt_values)/(float(alt_values) + float(ref_values))
						sample_data[info[0] + '*' +snp][sample][alt_codon] = frequency_alt
					else:
						sample_data[info[0] + '*' +snp][sample][alt_codon] = 1






# We read the .vcf files present in the folder that we introduce as an argument
if not args.path_to_vcf.endswith('/'):
	args.path_to_vcf += '/'
for vcf_file in os.listdir(args.path_to_vcf):
	if vcf_file.endswith('.vcf'):
		read_vcf(args.path_to_vcf, genome_positions, multistate_alleles, vcf_file, sample_data, samples, good_ones)



### Genome-wise pN/pS
all_syn = 0
all_non_syn = 0

for key in genome_positions.keys():
	for snp in genome_positions[key]:
		if genome_positions[key][snp].has_key('state2') and not multistate_alleles.has_key(key + '*' + snp):
			alt_codon = genome_positions[key][snp]['state1']
			ref_codon = genome_positions[key][snp]['state2']
			if gencode[alt_codon] == gencode[ref_codon]:
				all_syn += 1
			else:
				all_non_syn += 1
print('#pnps entire population: %s\n' %((float(all_non_syn)/float(genome_wise_possible['non_syn']))/(float(all_syn)/float(genome_wise_possible['syn']))))




### Per-gene table generation
output = open('%s_pnps_pergene.txt' %args.out_name, 'w')
output.write('gene\tsyn_found\tsyn_possible\tnon_syn_found\tnon_syn_possible\tgene_length\n')
for gene in genome_positions:
        non_syn_found = 0
        syn_found = 0
        for snp in genome_positions[gene]:
                if genome_positions[gene][snp].has_key('state2') and not multistate_alleles.has_key(gene + '*' + snp):
                        if gencode[genome_positions[gene][snp]['state1']] == gencode[genome_positions[gene][snp]['state2']]:
                                syn_found +=1
                        else:
                                non_syn_found += 1
                        reference = genome_positions[gene][snp]['state1']
        output.write('%s\t%s\t%s\t%s\t%s\t%s\n'%(gene, str(syn_found), str(per_gene_possible[gene]['syn']), str(non_syn_found), str(per_gene_possible[gene]['non_syn']),str(len(genome[gene]))))
output.close()
