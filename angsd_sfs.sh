#!/bin/sh

set -euxo pipefail

## angsd to calculate minor allele frequency (MAF) spectrum

## may require a lot of RAM!!! (up to 100 Gb)  

## usage: bash angsd_sfs.sh list_bam_files(full_path)($1)

## required software: angsd & scripts - realSFS & thetaStat from angsd package



	./angsd -bam ${1} -doSaf 1 -out output -anc reference_genome_23408contigs.fa -GL 2 -fold 1 -minMapQ 1 &&

	./realSFS output.saf.idx -maxIter 100 -P 10 > output.sfs


# OPTIONAL: calculates thetas, Tajima's D & other PopGen

# 	./angsd -bam ${1} -out output -doThetas 1 -doSaf 1 -anc reference_genome_23408contigs.fa -pest output.sfs -GL 2 -fold 1

#	./misc/thetaStat make_bed output.thetas.gz

# 	a=`gunzip -c output.thetas.gz | cut -f1 | sort -u | wc -l`

# 	./misc/thetaStat do_stat output.thetas.gz -nChr ${a}

# calculates average Tajima's D

# 	tail -n+2 output.thetas.gz.pestPG | awk '{sum +=$9;n++}END{print sum /n;}' > tajima.out
