# *Bathymodiolus Brooksi* symbionts

Welcome to BathyBrooksiSymbionts repository! Here, you will find the in-house scripts that were used for the [population structure and strain evolution](https://pages.github.com/) analyses of chemosynthetic symbionts from the deep-sea mussel *Bathymodiolus brooksi*. 

## Population_structure_analyses

In this folder are located the R and Python scripts that were used to study symbiont population structure. Additionally, you can find a toy example to test the scripts, which is composed by SNVs called in 6 different samples (A-F). 

### Simple example

First, we merge all SNPs found across samples into count files for each individual sample. This script is using .vcf files that have been produced by [**Lofreq**](http://csb5.github.io/lofreq/), and uses the DP4 field to extract the nucleotide counts.

`./vcf_to_mergedcounts.py .`
