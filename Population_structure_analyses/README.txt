Usage: Rscript structure.r <samples> <gff>

Example:

change into the folder data and run:
Rscript ../structure.r samples.txt Clue_final.gff

Input:

The samples file has one line for each sample with 2 fields: the sample name and filename of the corresponding SNP file.
Both the samples file and the SNP file must be in the current directory. however the R script can be in another directory and called with a path.
The SNP file must contain 7 fields: CHROM, Pos, Ref, A, C, G, and T.
The SNPs and the order of the SNPs MUST be the same in all files, the script does not check for this!

The gff file is a standard gff3 file.

Output:

<sample>Fst_pos.txt
<sample>Fst.txt

The Fst_pos.txt file contains one line for each SNP with the corresponding information information on the nucleotide counts in each sample and the sample and between-sample pi values.
The Fst.txt file contains one line for each gene with the corresponding pi values, Fst for the whole set of samples and pairwise Fst.
