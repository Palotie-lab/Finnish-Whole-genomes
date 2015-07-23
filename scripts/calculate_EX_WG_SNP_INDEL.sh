#!/bin/bash
### PROGRAM TO SUBSET .VCF FILES IN SNP/INDELS EXOME/NON EXOMES AND OUTPUT STATISTICS ###
### Argument 1: Input .vcf.gz file
### Argument 2: Output .vcf file subset, it also output two additional files *_gatkreport_* and *_bcfreport_* with stats
### Argument 3: It is "INDEL" if interested in indels or "SNP" if want to extract SNPs
### Argument 4: It is "YES" if want to extract only exomes and "NO" otherwise


source /home/unix/aganna/.my.bashrc
source /broad/software/scripts/useuse
use .zlib-1.2.8

#use Bcftools
#use Tabix

inD="$1"
outD="$2"
selectD="$3"
exomeD="$4"

if [ $exomeD == "YES" ]
    then 
        java -Xmx72g -jar /home/unix/aganna/GenomeAnalysisTK.jar \
        -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \
        -T SelectVariants \
        -V $inD \
        -selectType $selectD \
        -L /humgen/gsa-hpprojects/GATK/bundle/current/b37/Broad.human.exome.b37.interval_list \
        -o $outD
      
else 
    if [ $exomeD == "NO" ]
        then 
            java -Xmx180g -jar /home/unix/aganna/GenomeAnalysisTK.jar \
            -R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta \
            -T SelectVariants \
            -V $inD \
            -selectType $selectD \
            -o $outD
    else
        echo "The exome field should contain YES or NO"
    fi
fi

bgzip $outD
tabix -f $outD.gz