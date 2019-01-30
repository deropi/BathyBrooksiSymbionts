import sys

file1 = sys.argv[1]
file2 = sys.argv[2]

genes1 = set()
genes2 = set()

def read(file, genes):
	with open(file, 'r') as new_file:
		for line in new_file:
			genes.add(line.strip())
read(file1, genes1)
read(file2, genes2)

print len(genes1-genes2)
