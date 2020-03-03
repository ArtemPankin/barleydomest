## Scripts accompanying New Phytologist 2018 manusctipt

**MODULE 1 (2 files) extracting putative ancestry of the domesticated genes using the Maxumum Likelihood distance approach**

1. run gene_ancestry_raxml.sh bash script
2. use output of the bash script in gene_ancestry.R

**MODULE 2 (3 files): read filtering, mapping, SNP calling and filtering; parallelized pipeline (bsub LSF)**

Starting script:
read_mapping_snp_calling.sh

Companion scripts should be in the folder, bsub LSF system required:
filt_map_snp.bsub
genotyping.bsub


**MODULE 3: calculating site frequency spectra (sfs) - a.k.a. minor allele frequency (MAF) - using two methods - based on SNP genotypes in a vcf file (1) and on raw bam files using angsd (2)**

sfs_angsd_vs_vcf.R - R script to plot and compare sfs spectra

**ADDITIONAL FILES**

reference_genome_23408contigs.fa - reference genome for the read mapping

cDNA_coordinates.bed - list of the reference targets present in the cDNA form (required for trimmming of the gDNA reads that are mapped onto a cDNA reference)

wild_dom_no_adm90 - annotation file for the non-admixed wild and domesticated genotypes 
