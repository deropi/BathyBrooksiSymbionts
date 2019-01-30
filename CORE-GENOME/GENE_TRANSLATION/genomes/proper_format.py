import sys
fasta = sys.argv[1]
with open(fasta, 'r') as input_file, open('formated_%s'%fasta, 'w') as output_file:
	for line in input_file:
		if line.startswith('>'):
			header = line.strip()
			seq = input_file.next().strip()
			output_file.write(header + '\n')
			for i in range(0, len(seq), 60):
				output_file.write(seq[i:i+60] + '\n')
