#!/bin/bash
### PROGRAM TO CREATE CONCORDANCE STATS BETWEEN WGS DATA AND CHIP DATA  ###
### Argument 1: Input file in plink (e.g. G77318RH) containing chip data.
### Argument 2: Input .vcf.gz file from WGS to be compared.
### Argument 3: Name of the output analysis (from GATK)

## imput variables ##
inchip=$1
inWGS=$2
outD=$3

#inchip="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/FinnRisk"
#inWGS="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_SNP_EX.vcf.gz"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_EX"

## Save positions and SNP id name for all the chip SNPs, excluding problematic SNPs ##
awk '$4 != "-" && $5 != "-" && $1 != "0" && $2 !~ /^indel/ && $2 !~ /^1KG/ {print $1,$4,$4,$2 }'  ${inchip}.bim  | sed 's/ /\t/g' > /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp1

## Save position for all the variants in the WGS file ##
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\n" $inWGS > /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp2

## Merge the two files based on CHR and positions, keep only positions in common ##
## Sort file because unsorted my create some problem in the subsetting of the files ##

awk 'NR==FNR{a[$1,$2]=$3 FS $4;next} ($1,$2) in a{print $0, a[$1,$2]}' /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp2 /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp1 | sort -n --key=1,1 --key=2,2  > /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3


## ASSIGN CORRECT ALLELES ##
/home/unix/aganna/plink_linux_x86_64/plink -bfile $inchip \
--extract range /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3 \
--a2-allele /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3 5 4 \
--snps-only \
--list-duplicate-vars ids-only suppress-first \
--real-ref-alleles \
--make-bed \
--out /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp4

## Remove again indels which might have pop-up ##
cat  /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp4.dupvar > /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3ind
awk '$2 ~ /^indel/ {print $2 }'  /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp4.bim  >> /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3ind

## And now recode to .vcf ##
/home/unix/aganna/plink_linux_x86_64/plink -bfile /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp4 \
--exclude /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3ind \
--a1-allele /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3 6 4 \
--geno 0.05 \
--recode vcf-iid \
--out /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp5


## Zip data ##
bgzip -f /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp5.vcf
tabix /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp5.vcf.gz


## Extract regions from the WG file
bcftools norm -D -R  /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp3 $inWGS -Oz -o /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.tempFWG.vcf.gz

tabix -f /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.tempFWG.vcf.gz

java -Xmx68g -jar /home/unix/aganna/GenomeAnalysisTK.jar \
     -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \
     -T GenotypeConcordance \
     --moltenize \
     -eval /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.tempFWG.vcf.gz \
     -comp /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp5.vcf.gz \
     -o ${outD}_DP_ALL.txt
         
for i in {5..40..5}
  do
     eval "java -Xmx68g -jar /home/unix/aganna/GenomeAnalysisTK.jar \
     -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \
     -T GenotypeConcordance \
     --moltenize \
     -eval /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.tempFWG.vcf.gz \
     -comp /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp5.vcf.gz \
     -gfe 'DP > $i' \
     -o ${outD}_DP_${i}.txt"
 done

## Remove temporary files ##
rm /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/.temp*

