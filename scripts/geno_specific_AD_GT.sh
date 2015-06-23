#!/bin/bash
### PROGRAM TO OBTAIN GENOTYPE-SPECIFIC AD. IT IS MEANINGFULL ONLY FOR SNPs ###
### Argument 1: Input .vcf or .vcf.gz file
### Argument 2: Output .txt with the results


source /home/unix/aganna/.my.bashrc

inD="$1"
outD="$2"

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_SNP_EX.vcf.gz"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/AD_G77318RH_SNP_EX.txt"


bcftools view -g het -v snps $inD |
  sed 's/^##FORMAT=<ID=AD,Number=\./##FORMAT=<ID=AD,Number=R/g' |
  bcftools norm -Ou -m -any |
  bcftools norm -Ou -f /humgen/1kg/reference/human_g1k_v37.fasta | bcftools query -f "[%GT %AD ]\n" > $outD

   

  