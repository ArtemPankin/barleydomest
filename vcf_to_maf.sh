#!/bin/sh

# This scripts calculates numeric stats for each SNP

# output 10 columns: conitg	position	missing_count	reference_count	heterozygous_count	homozygous_count	heterozygous_fraction	missing_fraction	hom_het_fraction	homozygous_fraction

# usage: bash vcf_to_maf.sh input.vcf($1)

grep -v "#" $1 | 
awk '
{
	for(i=1;i<=NF;i++)
		if($i ~ /^1\/1\:/)hom++; 
		else if($i ~ /^0\/1\:/)het++; 
		else if($i ~ /^0\/0\:/)ref++; 
		else if($i ~ /^.\/\./)mis++
	}
	{
		print $1 "\t" $2 "\t" mis "\t" ref "\t" het "\t" hom "\t" het/(hom+het+ref) "\t" mis/(hom+het+ref+mis) "\t" (hom*2+het)/(2*(het+hom+ref)) "\t" hom/(hom+ref); hom=0; het=0; ref=0; mis=0
	}' | 
awk '
BEGIN{
	FS = OFS = "\t"
	}
	{
		for(i=1;i<=NF;i++)
			if($i ~ /^$/)$i = 0
			};1' | 
awk '
{
	if($NF>0.5)print$0"\t"1-$NF; 
	else print$0"\t"$NF
	}' > $1.snp_stats.out
