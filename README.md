## Maria von Korff group (MPIPZ / HHU / CEPLAS) script repository


### MODULE (2 files) extracting putative ancestry of the domesticated genes using the Maxumum Likelihood distance approach

1. run gene_ancestry_raxml.sh bash script
2. use output of the bash script in gene_ancestry.R

Generates figures Fig. 2de

### MODULE (3 files): read filtering, mapping & SNP calling and filtering

Starting script:
read_mapping_snp_calling.sh

Companion scripts should be in the folder, bsub LSF system required:
filt_map_snp.bsub
genotyping.bsub

### ADDITIONAL FILES

reference_genome_23408contigs.fa - reference genome for the read mapping

cDNA_coordinates.bed - list of the reference targets present in the cDNA form (required for trimmming of the gDNA reads that are mapped onto a cDNA reference)

wild_dom_no_adm90 - annotation file for the non-admixed wild and domesticated genotypes 
