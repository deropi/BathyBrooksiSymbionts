# *Bathymodiolus Brooksi* symbionts

Welcome to BathyBrooksiSymbionts repository! Here, you will find the in-house scripts that were used for the [population structure and strain evolution](https://pages.github.com/) analyses of chemosynthetic symbionts from the deep-sea mussel *Bathymodiolus brooksi*. 

## Population structure analysis

This folder contains R and Python scripts that were used to study symbiont population structure and selective preassure. Additionally, you can find a toy example to test the scripts, which is composed by SNVs called in 6 different samples (A-F). Note that these analyses are gene-based, thus, the reference sequences and coordinates for the SNVs are based on genes.

### Simple example

First, we merge all SNPs found across samples into count files for each individual sample. This script is using .vcf files that have been produced by [**Lofreq**](http://csb5.github.io/lofreq/), and requires the DP4 field to extract the nucleotide counts.

`./vcf_to_mergedcounts.py .`

This will produce `.list` files for each of the samples. A two columns `.txt` file is required, where the name of the samples must be linked to the name of its corresponding file (see `metagenomes/samples.txt`). This table is used as the input for the R script, which will perform the population structure analyses. 

`Rscript structure.r samples.txt`

This results in the generation of a matrix `samples.txtFst_pos.txt` ,containing Fst and Pi values for each SNV found in each of the samples. Then, we use this table as input to `genome_wise_calculations.py`, which will estimate Fst and Pi genome-wide.

`./genome_wise_calculations.py samples.txtFst_pos.txt reference_genome.fasta toy_example`

This will output a dataframe with Pi values for each sample (`toy_example_pi.out`), and a matrix with pair-wise Fst values (`toy_example_fst.mat`). 

Additionally, pN/pS per gene and for the entire population are estimated with the following command:

`./pN_pS_calculations.py . reference_genome.fasta toy_pNpS good_snps.txt`

A file with the SNVs that must be parsed is required. This file contains in each line a SNV ID, which contains the geneID and the position in the following format: `Gene1*3`.

