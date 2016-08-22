#!/bin/sh

## This scripts calculates pairwise ML distances between each reference targets (1st column of vcf), genes and contigs in this case
## Sofware in the bin path: raxmlHPC-SSE3
## File 'wild_dom_no_adm90' (available on github) should be in the folder
## It generates an input file for the 'gene_ancestry.R' script
## Usage: bash gene ancestry_raxml.sh input.vcf


# extract fasta files for each gene > 20 SNPs

vcf=$1

	mkdir -p {gene_by_gene,raxml_dist}

## set number of samples variable

samnum=`tail -n 1 ${vcf} | awk '{print NF}'`

## extract list of targets/genes with > 20 SNPs 

	grep -v "#" ${vcf} | cut -f 1 | uniq -c | head -n 100 | awk '$1 > 20' > gene.list.20

## convert vcf to fa for each gene

while read line; do ## while1 START 

if [[ ! -s gene_by_gene/$line.$vcf.fa ]]; then # if1 START

	grep -w "$line" ${vcf} > temp.$line.$vcf 

for ((i=10;i<=samnum;i++)); do # for1 START

	cat <(echo ">"`cut -f${i} ids.noout`"_"${line}) \
	<(echo `cat temp.$line.$vcf | 
	cut -f4,5,${i} | 
	tail -n+2 | 
	awk '{if($3 ~ "0/0")print $1; else if($3 ~ "1/1")print $2; else if($3 ~ "./.")print "-"}' | 
	tr -d "\n"`)

done > gene_by_gene/$line.$vcf.fa # for1 END

fi # if1 END

	rm -f temp.$line.$vcf

##

## calculate ML distances 

if [[ ! -s raxml_dist/RAxML_distances.${line}.temp ]]; then # if2 START

	./raxmlHPC-SSE3 -f x -p 12345 -s gene_by_gene/$line.$vcf.fa -m GTRGAMMA -n $line.temp

	mv RAxML_distances.${line}.temp raxml_dist/

	rm -f RAxML_parsimonyTree.${line}.temp 

	mv RAxML_info.${line}.temp raxml_dist/

fi # if2 END

## select results only for wild vs. domesticated genotypes and annotate with genotype name and wild/domesticated status

if [[ ! -s gene_by_gene/distances/dist.$line ]]; then #if3 START

	awk 'NR==FNR{a[$1]=$10;next}{if($1 in a)print a[$1]"\t"a[$2]"\t"$0}' wild_dom_no_adm90 \
	<(cat raxml_dist/RAxML_distances.${line}.temp | 
	sed "s/_${line}//g") | 
	egrep -v "wild.*wild|dom.*dom" | 
	awk 'BEGIN{OFS="\t"}{if($1 ~ "wild")print $2,$1,$4,$3,$5,"'${line}'"; else print $0,"'${line}'"}' | sort -k3,3 -k5,5n > raxml_dist/dist.$line
	
fi # if3 END

	
	cat raxml_dist/* > MLdist.data


done < gene.list.20 # while1 END
