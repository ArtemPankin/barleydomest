#!/bin/sh

set -euxo pipefail

## This script performs Illumina read filtering, mapping and SNP calling and filtering

## The following SOFTWARE should be in the 'bin' path: samtools, bwa, GATK, cd-hit-dup, lighter, intersectBed, bamToFastq, AddOrReplaceReadGroups.jar, 
## The script requires two dependent BSUB scripts 'filt_map_snp.bsub', 'genotyping.bsub'. They should be placed in the same folder
## we recommend reserving memory for this script (~ 80 Gb in out case)

## The script should be started in with the eight compulsory arguments: 
## bash read_mapping_snp_calling.sh /
## $1 - reference_file (.fa) / 
## $2 - list_of_fastq files (one file name per line, no extra characters/FULL_PATH) /
## $3 - number of cores for multicore processing /
## $4 - dataset type filtered/unfiltered ('filtered' - omits filtering steps; 'unfiltered' - proceeds with filtering) /
## $5 - coordinates of the cDNA reference in bed format /
## $6 - first base to keep in bp, based on the FASTQC results / 
## $7 - last base to keep in bp, based on the FASTQC results /
## $8 - error correction by Lighter software (correction/nocorrection) /
## $9 - skip_mapping/skip_genotyping(optional)

## Example command:
## bash read_mapping_snp_calling.sh reference_genome_23408contigs.fa read_files.list 8 unfiltered cDNA_coordinates.bed 5 90 correction

## The script is using BSUB LSF system for parallel computing

#!/bin/sh

# counts number of fastq files to process

if [ $4 = "unfiltered" ]; then

	grep -v "_2.fq" $2 > _list

	num=`awk '{sub("\r$","");print}' _list | wc -l`

	
# modifies the .bsub script

	sed "s/1-.*]/1-$num]/" filt_map_snp.bsub | sed -r "s/BSUB -n .*$/BSUB -n $3/" | sed "s:\$1:$1:g" | sed "s:\$2:_list:g" | sed "s:\$3:$3:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$6:$6:g" | sed "s:\$7:$7:g" | sed "s:\$8:$8:g" | awk '{sub("\r$","");print}' > filt_map_snp.temp.bsub

else

	num=`awk '{sub("\r$","");print}' $2 | wc -l`

	sed "s/1-.*]/1-$num]/" filt_map_snp.bsub | sed -r "s/BSUB -n .*$/BSUB -n $3/" | sed "s:\$1:$1:g" | sed "s:\$2:$2:g" | sed "s:\$3:$3:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$6:$6:g" | sed "s:\$7:$7:g" | sed "s:\$8:$8:g" | awk '{sub("\r$","");print}' > filt_map_snp.temp.bsub

fi

# indexing reference if the index files do not exist, else skip

if [ ! -f ${1}.fai ]; then 

	samtools faidx $1 

fi

if [ ! -f ${1}.bwt ]; then 

	bwa index $1

fi

if [ ! -s ${1%%.fa}.dict ]; then

	java -jar CreateSequenceDictionary.jar R= $1 O= ${1%%fa}dict

fi

if [ ! ${9} = "skip_mapping" ]; then


bsub < filt_map_snp.temp.bsub > process_id


# sleeping until all bjobs are finished; comparing DONE and JOBS fields in 'bjobs -A' output

	subm=`awk -F "<|>" '{print$2}' process_id | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished=`awk -F "<|>" '{print$2}' process_id | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited=`awk -F "<|>" '{print$2}' process_id | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed=$[$finished+$exited]

while [ $subm != $completed ]; do

	sleep 20s 

	subm=`awk -F "<|>" '{print$2}' process_id | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished=`awk -F "<|>" '{print$2}' process_id | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited=`awk -F "<|>" '{print$2}' process_id | grep -f - <(bjobs -A) | awk '{print$8}'`
	
	completed=$[$finished+$exited]

	echo $subm
	echo $completed

done

	rm process_id


# verificiation of SNP calling step

	grep -f <(ls -l *raw.vcf.idx| grep -w "0" | awk -F '_| ' '{print$(NF-1)"_"}') $2 > temp.failed.list

while [ -s temp.failed.list ]; do

if [ $4 = "unfiltered" ]; then

	grep -v "_2.fq" temp.failed.list > temp_failed_sample_list

	num1=`awk '{sub("\r$","");print}' temp_failed_sample_list | wc -l`

# modifies the .bsub script

	sed "s/1-.*]/1-$num1]/" filt_map_snp.bsub | sed -r "s/BSUB -n .*$/BSUB -n $3/" | sed "s:\$1:$1:g" | sed "s:\$2:temp_failed_sample_list:g" | sed "s:\$3:$3:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$6:$6:g" | sed "s:\$7:$7:g" | sed "s:\$8:$8:g" | awk '{sub("\r$","");print}' > filt_map_snp.temp1.bsub

else

	num1=`awk '{sub("\r$","");print}' temp.failed.list | wc -l`

	sed "s/1-.*]/1-$num1]/" filt_map_snp.bsub | sed -r "s/BSUB -n .*$/BSUB -n $3/" | sed "s:\$1:$1:g" | sed "s:\$2:temp.failed.list:g" | sed "s:\$3:$3:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$6:$6:g" | sed "s:\$7:$7:g" | sed "s:\$8:$8:g" | awk '{sub("\r$","");print}' > filt_map_snp.temp1.bsub

fi

# submitting failed samples SNP calling bsub script

	bsub < filt_map_snp.temp1.bsub > process_id2

# sleeping until all bjobs are finished; comparing DONE and JOBS fields in 'bjobs -A' output

	subm3=`awk -F "<|>" '{print$2}' process_id2 | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished3=`awk -F "<|>" '{print$2}' process_id2 | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited3=`awk -F "<|>" '{print$2}' process_id2 | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed3=$[$finished3+$exited3]

while [ $subm3 != $completed3 ]

do

	sleep 20s 

	subm3=`awk -F "<|>" '{print$2}' process_id2 | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished3=`awk -F "<|>" '{print$2}' process_id2 | grep -f - <(bjobs -A) | awk '{print$6}'`
	
	exited3=`awk -F "<|>" '{print$2}' process_id2 | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed3=$[$finished3+$exited3]
	

	echo $subm3

	echo $completed3

done

	rm process_id2

grep -f <(ls -l *raw.vcf.idx| grep -w "0" | awk -F '_| ' '{print$(NF-1)"_"}') $2 > temp.failed.list

done

fi ## closing "skip_mapping"

# merging multiple VCFs to get the overall set of unfiltered SNPs and INDELs

	raw_snps=`ls -l *raw.vcf | grep -f <(awk '{FS="/|_"}{print $(NF-1)}' $2) | awk '{print$NF}' | awk 'BEGIN{RS=""}{gsub("\n"," ",$0);print}' | awk 'BEGIN{RS=""}{gsub("\n"," ",$0);print}' | awk '{gsub(" "," --variant ");print "--variant " $0}'` 


if [ ! -s rawcalls.${num}samples.$2.UG.minimumN1.vcf.idx ]; then

	java -Xmx80g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out rawcalls.${num}samples.$2.UG.minimumN1.vcf $raw_snps --genotypemergeoption UNSORTED

fi

## opening "skip_genotyping" condition

if [ ${9} = "skip_genotyping" ] || [ ${10} = "skip_genotyping" ]; then

break

else

# submitting genotyping bsub script

if [ $4 = "unfiltered" ]; then


	grep -v "_2.fq" $2 > _list

	num=`awk '{sub("\r$","");print}' _list | wc -l`

# modifies the .bsub script

	sed "s:\$1:$1:g" genotyping.bsub | sed "s:1-.*]:1-$num]:" | sed "s:\$2:_list:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$8:$8:g" | sed "s:\$9:rawcalls.${num}samples\.$2\.UG\.minimumN1\.vcf:g" | sed "s:\$num:$num:g" > genotyping.temp.bsub

else

	num=`awk '{sub("\r$","");print}' $2 | wc -l`

	sed "s:\$1:$1:g" genotyping.bsub | sed "s:1-.*]:1-$num]:" | sed "s:\$num:$num:g" | sed "s:\$2:$2:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$8:$8:g" | sed "s:\$9:rawcalls.${num}samples\.$2\.UG\.minimumN1\.vcf:g" > genotyping.temp.bsub


fi

	bsub < genotyping.temp.bsub > process_id1

# sleeping until all bjobs are finished; comparing DONE and JOBS fields in 'bjobs -A' output

	subm1=`awk -F "<|>" '{print$2}' process_id1 | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished1=`awk -F "<|>" '{print$2}' process_id1 | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited1=`awk -F "<|>" '{print$2}' process_id1 | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed1=$[$finished1+$exited1]

while [ $subm1 != $completed1 ]

do

	sleep 20s 

	subm1=`awk -F "<|>" '{print$2}' process_id1 | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished1=`awk -F "<|>" '{print$2}' process_id1 | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited1=`awk -F "<|>" '{print$2}' process_id1 | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed1=$[$finished1+$exited1]
	
	echo $subm1

	echo $completed1

done

	rm process_id1



# verificiation of SNP genotyping step
	
	grep -f <(ls -l *SNP.homo.UG.genotyping.DP8.pcrerr-5e-2.${num}.vcf | grep -w "0" | awk -F '_| ' '{print$(NF-1)"_"}') $2 > temp1.failed.list

while [ -s temp1.failed.list ]; do

if [ $4 = "unfiltered" ]; then

	grep -v "_2.fq" temp1.failed.list > temp1_failed_sample_list

	num2=`awk '{sub("\r$","");print}' temp1_failed_sample_list | wc -l`

# modifies the .bsub script

	sed "s:\$1:$1:g" genotyping.bsub | sed "s:1-.*]:1-$num2]:" | sed "s:\$2:temp1_failed_sample_list:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$8:$8:g" | sed "s:\$9:rawcalls.${num}samples\.$2\.UG\.minimumN1\.vcf:g" | sed "s:\$num:$num:g" > genotyping.temp1.bsub

else

	num2=`awk '{sub("\r$","");print}' temp1.failed.list | wc -l`

	sed "s:\$1:$1:g" genotyping.bsub | sed "s:1-.*]:1-$num2]:" | sed "s:\$2:temp1.failed.list:g" | sed "s:\$4:$4:g" | sed "s:\$5:$5:g" | sed "s:\$8:$8:g" | sed "s:\$9:rawcalls.${num}samples\.$2\.UG\.minimumN1\.vcf:g" | sed "s:\$num:$num:g" > genotyping.temp1.bsub

fi

# submitting failed samples SNP calling bsub script

	bsub < genotyping.temp1.bsub > process_id3

# sleeping until all bjobs are finished; comparing DONE and JOBS fields in 'bjobs -A' output

	subm2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed2=$[$finished2+$exited2]

while [ $subm2 != $completed2 ]

do

	sleep 20s 

	subm2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$4}'`

	finished2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$6}'`

	exited2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$8}'`

	exited2=`awk -F "<|>" '{print$2}' process_id3 | grep -f - <(bjobs -A) | awk '{print$8}'`

	completed2=$[$finished2+$exited2]

	echo $subm2

	echo $completed2

done

	rm process_id3

	grep -f <(ls -l *SNP.homo.UG.genotyping.DP8.pcrerr-5e-2.${num}.vcf | grep -w "0" | awk -F '_| ' '{print$(NF-1)"_"}') $2 > temp1.failed.list

done


fi ## closing "skip_genotyping" condition

# when finished merge genotyping into multisample vcf

	genot_snps_hom_het=`ls -l *SNP.homo_het*pcrerr-5e-2.$num.vcf | grep -f <(awk '{FS="/|_"}{print $(NF-1)}' $2) | awk '{print$NF}' | awk 'BEGIN{RS=""}{gsub("\n"," --variant ");print "--variant " $0}'`

	genot_indels_hom_het=`ls -l *INDEL.homo_het*pcrerr-5e-2.$num.vcf | grep -f <(awk '{FS="/|_"}{print $(NF-1)}' $2) | awk '{print$NF}' | awk 'BEGIN{RS=""}{gsub("\n"," --variant ");print "--variant " $0}'`

	genot_snps_hom=`ls -l *SNP.homo.*pcrerr-5e-2.$num.vcf | grep -f <(awk '{FS="/|_"}{print $(NF-1)}' $2) | awk '{print$NF}' | awk 'BEGIN{RS=""}{gsub("\n"," --variant ");print "--variant " $0}'`

	genot_indels_hom=`ls -l *INDEL.homo.*pcrerr-5e-2.$num.vcf | grep -f <(awk '{FS="/|_"}{print $(NF-1)}' $2) | awk '{print$NF}' | awk 'BEGIN{RS=""}{gsub("\n"," --variant ");print "--variant " $0}'`


# merging all VCFs to get the overall set of filtered SNPs and INDELs

if [ ! -s $2.SNP.homo.${num}samples.UG.genotyping.minN1.raw.vcf.idx ]; then

	java -Xmx80g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.SNP.homo.${num}samples.UG.genotyping.minN1.raw.vcf $genot_snps_hom --genotypemergeoption UNSORTED

fi

if [ ! -s $2.INDEL.homo.${num}samples.UG.genotyping.minN1.raw.vcf.idx ]; then

	java -Xmx50g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.INDEL.homo.${num}samples.UG.genotyping.minN1.raw.vcf $genot_indels_hom --genotypemergeoption UNSORTED

fi

if [ ! -s $2.SNP.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf.idx ]; then

	java -Xmx80g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.SNP.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf $genot_snps_hom_het --genotypemergeoption UNSORTED

fi

if [ ! -s $2.INDEL.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf.idx ]; then

	java -Xmx50g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.INDEL.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf $genot_indels_hom_het --genotypemergeoption UNSORTED

fi

## various filtering steps, optional

:<<"COMMENT"

# filtering out non-variant positions and REMOVING singleton SNPs and INDELS after merging the filtered files

# hom_het

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.wo_sing.vcf.idx ]; then

	java -Xmx10g -jar GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.SNP.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 2' -o $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.wo_sing.vcf

fi

if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.wo_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.INDEL.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 2' -o $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.wo_sing.vcf

fi

# only hom

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.wo_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.SNP.homo.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 2' -o $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.wo_sing.vcf

fi

if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.wo_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.INDEL.homo.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 2' -o $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.wo_sing.vcf

fi

# filtering out non-variant positions and LEAVING singleton SNPs and INDELS after merging the filtered files

# hom_het

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.with_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.SNP.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 1' -o $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.with_sing.vcf

fi

if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.with_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.INDEL.homo_het.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 1' -o $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.with_sing.vcf

fi

# only hom

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.SNP.homo.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 1' -o $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.vcf

fi

if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.with_sing.vcf.idx ]; then

	java -Xmx10g -jar /biodata/dep_coupland/grp_korff/bin/GenomeAnalysisTK.jar -R $1 -T SelectVariants --variant $2.INDEL.homo.${num}samples.UG.genotyping.minN1.raw.vcf -select 'AC > 1' -o $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.with_sing.vcf

fi


# selecting sites with at least 1 homozygous SNP (with and wo sing); only for homo_het

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.wo_sing.hetless100.vcf ]; then

	grep -E "#|1\/1" $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.wo_sing.vcf > $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.wo_sing.hetless100.vcf

fi

if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.wo_sing.hetless100.vcf ]; then

	grep -E "#|1\/1" $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.wo_sing.vcf > $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.wo_sing.hetless100.vcf

fi

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.with_sing.hetless100.vcf ]; then

	grep -E "#|1\/1" $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.with_sing.vcf > $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.with_sing.hetless100.vcf	

fi

if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.with_sing.hetless100.vcf ]; then

	grep -E "#|1\/1" $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.with_sing.vcf > $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.with_sing.hetless100.vcf

fi

# selecting HOMOZYGOUS SNPs without missing data

# with sing

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.with_sing.vcf ]; then

grep -v "\.\/\." $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.vcf > $2.filtered.SNP.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.with_sing.vcf 	

fi

if [ ! -s $2.filtered.no_missing.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.with_sing.vcf ]; then

grep -v "\.\/\." $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.with_sing.vcf > $2.filtered.INDEL.${num}samples.UG.genotyping.no_missing.DP10.pcrerr5e-2.homo.with_sing.vcf

fi

# wo sing

if [ ! -s $2.filtered.SNP.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.wo_sing.vcf ]; then

grep -v "\.\/\." $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.wo_sing.vcf > $2.filtered.SNP.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.wo_sing.vcf 	

fi


if [ ! -s $2.filtered.INDEL.${num}samples.UG.genotyping.no_missing.DP10.pcrerr5e-2.homo.wo_sing.vcf ]; then

grep -v "\.\/\." $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo.wo_sing.vcf > $2.filtered.INDEL.${num}samples.UG.genotyping.no_missing.DP10.pcrerr5e-2.homo.wo_sing.vcf

fi

# merging SNPs and INDELs

# hom + het

	java -Xmx20g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.fully_filtered.genotyping.${num}samples.INDEL.SNP.homo_het.wo_sing.vcf --variant $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.wo_sing.hetless100.vcf --variant $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.wo_sing.hetless100.vcf --genotypemergeoption UNSORTED

	java -Xmx20g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.fully_filtered.genotyping.${num}samples.INDEL.SNP.homo_het.with_sing.vcf --variant $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo_het.with_sing.hetless100.vcf --variant $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP10.pcrerr5e-2.homo_het.with_sing.hetless100.vcf --genotypemergeoption UNSORTED

# only hom

	java -Xmx20g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.fully_filtered.genotyping.${num}samples.INDEL.SNP.homo.wo_sing.vcf --variant $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.wo_sing.vcf --variant $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.wo_sing.vcf --genotypemergeoption UNSORTED

	java -Xmx20g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.fully_filtered.genotyping.${num}samples.INDEL.SNP.homo.with_sing.vcf --variant $2.filtered.INDEL.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.vcf --variant $2.filtered.SNP.${num}samples.UG.genotyping.minN1.DP8.pcrerr5e-2.homo.with_sing.vcf --genotypemergeoption UNSORTED

# only hom + no_missing

	java -Xmx20g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.fully_filtered.no_missing.genotyping.${num}samples.INDEL.SNP.homo.wo_sing.vcf --variant $2.filtered.SNP.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.wo_sing.vcf --variant $2.filtered.INDEL.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.wo_sing.vcf --genotypemergeoption UNSORTED

	java -Xmx20g -jar GenomeAnalysisTK.jar -R $1 -T CombineVariants --minimumN 1 --out $2.fully_filtered.no_missing.genotyping.${num}samples.INDEL.SNP.homo.with_sing.vcf --variant $2.filtered.SNP.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.with_sing.vcf --variant $2.filtered.INDEL.${num}samples.UG.genotyping.no_missing.DP8.pcrerr5e-2.homo.with_sing.vcf --genotypemergeoption UNSORTED

COMMENT


rm _list
rm genotyping.temp.bsub
rm filt_map_snp.temp.bsub
rm temp.failed.list
rm temp_failed_sample_list
rm filt_map_snp.temp1.bsub
rm genotyping.temp.bsub
rm temp1.failed.list
rm temp1_failed_sample_list
