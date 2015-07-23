#!/bin/bash
### PROGRAM TO CALCULTE KINSHIP MATRIX, MDS AND SEX CHECK  ###
### Argument 1: Input file in plink (e.g. G77318RH). DO NOT specify the extension.
### Argument 2: Folder where to output the results

source /home/unix/aganna/.my.bashrc

inD="$1"
outD="$2"

## This script is used to obtained high-confidence SNPs from the 1000genome project ##
awk  '$1 ~ /^[^#]/ { print $1, " " , $2, " " , $2, " ",$3  }' /humgen/gsa-hpprojects/GATK/bundle/current/b37/1000G_phase1.snps.high_confidence.b37.vcf  > ${outD}positions_to_include.txt

## Add extract regions to exclude ##
echo -e "6\t25800000\t36000000\t.\n8\t6000000\t16000000\t.\n17\t40000000\t40000000\t." > ${outD}positions_to_exclude.txt

## Sex-check
#/home/unix/aganna/plink_linux_x86_64/plink -bfile ${inD} \
#--maf 0.05 \
#--extract range ${outD}positions_to_include.txt \
#--check-sex \
#--hwe  0.000001 \
#--geno 0.05 \
#--out ${outD}


## Exclude regions and clean the data ##

/home/unix/aganna/plink_linux_x86_64/plink -bfile $inD \
--extract range ${outD}positions_to_include.txt \
--make-bed \
--out ${outD}_QC1


/home/unix/aganna/plink_linux_x86_64/plink -bfile ${outD}_QC1 \
--exclude range ${outD}positions_to_exclude.txt \
--snps-only \
--hwe  0.000001 \
--maf 0.05 \
--geno 0.05 \
--make-bed \
--out ${outD}_QC2

## Prune the data ##
/home/unix/aganna/plink_linux_x86_64/plink -bfile ${outD}_QC2 \
--indep-pairwise 100 50 0.2 \
--out ${outD}_QC2

## Keep only pruned data ##
/home/unix/aganna/plink_linux_x86_64/plink -bfile ${outD}_QC2 \
--extract ${outD}_QC2.prune.in \
--make-bed \
--out ${outD}_QC2_pruned


## Run IBD ##
king -b ${outD}_QC2_pruned.bed --kinship --prefix ${outD}_QC2_pruned

## Run MDS ##
king -b ${outD}_QC2_pruned.bed --mds --prefix ${outD}_QC2_pruned


