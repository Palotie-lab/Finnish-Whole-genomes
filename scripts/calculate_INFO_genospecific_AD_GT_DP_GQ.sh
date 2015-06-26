#!/bin/bash
### PROGRAM TO OBTAIN TABLE-FORMAT RESULTS FOR EACH VARIANT AND EXPORT GENOTYPE SPECIFIC DP, GQ, GT and AD ###
### Argument 1: Input .vcf or .vcf.gz file
### Argument 2: folder where to save the results
### Argument 3: Output suffix for the four files. INFO, DP, GQ, ADGT will be added in front of the suffix. 
###             if nrow of the file is > 10M then 10M rows are sampled (true only for DP and GQ files)

source /home/unix/aganna/.my.bashrc

inD="$1"
outD="$2"
suffixD="$3"

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_SNP_EX.vcf.gz"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/"
#suffixD="G77318RH_SNP_EX"


java -Xmx72g -jar /home/unix/aganna/GenomeAnalysisTK.jar -T VariantEval \
-R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta  \
-o ${outD}COUNT_gatkreport_${suffixD}.txt \
--eval:set1 $inD \
-ST Sample -noST -noEV -EV CountVariants

java -Xmx72g -jar /home/unix/aganna/GenomeAnalysisTK.jar -T VariantEval \
-R /humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta  \
-o ${outD}TITV_gatkreport_${suffixD}.txt \
--eval:set1 $inD \
-ST Sample -noST -noEV -EV TiTvVariantEvaluator

bcftools stats $inD > ${outD}bcfreport_${suffixD}.txt


bcftools norm -Ou  -m -any $inD | bcftools query -Hf '%CHROM %POS %ID %REF %ALT %FILTER %AF %QD %DP %MQ %GQ_MEAN %InbreedingCoeff %FS %MQRankSum %ReadPosRankSum %AC %VQSLOD\n'  > ${outD}INFO_${suffixD}.txt

bcftools view -g het -v snps $inD |
  sed 's/^##FORMAT=<ID=AD,Number=\./##FORMAT=<ID=AD,Number=R/g' |
  bcftools norm -Ou -m -any |
  bcftools norm -Ou -f /humgen/1kg/reference/human_g1k_v37.fasta | bcftools query -f "[%GT %AD ]\n" > ${outD}ADGT_${suffixD}.txt
      
nlines=`wc -l $outD | awk '{print $1}'`

if [ $nlines -gt 10000000 ]; then
    echo 'Subsampling only 10M rows'
    bcftools norm -Ou  -m -any $inD | bcftools query -f '[%DP ]\n' | shuf -n 10000000 > ${outD}DP_${suffixD}.txt
    bcftools norm -Ou  -m -any $inD | bcftools query -f '[%GQ ]\n' | shuf -n 10000000 > ${outD}GQ_${suffixD}.txt
else
    bcftools norm -Ou  -m -any $inD | bcftools query -f '[%DP ]\n' > ${outD}DP_${suffixD}.txt
    bcftools norm -Ou  -m -any $inD | bcftools query -f '[%GQ ]\n' > ${outD}GQ_${suffixD}.txt
fi

    