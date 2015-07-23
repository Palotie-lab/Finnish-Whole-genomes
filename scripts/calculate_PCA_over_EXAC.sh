#!/bin/bash
### PROGRAM TO PCA AND OVERLAP WITH EXAC  ###
### Argument 1: Input file in vcf.gz format. 
### Argument 2: Folder where to output the results (.csv file containing the PCs)
### Notice that the script uses the /temp/ folder and all the scripts should be stored in /home/unix/aganna/scripts/PCA/ 

source /home/unix/aganna/.my.bashrc
source /broad/software/scripts/useuse
use .zlib-1.2.8

inD=$1
outD=$2

#inD=/humgen/atgu1/fs03/wip/aganna/BDS_exomes/original/seq/C1836_PASS_NLC.vcf.gz
#outD=/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/

/home/unix/aganna/scripts/PCA/extract_sites.pl ${inD} /broad/hptmp/aganna/outputPCA.vcf /home/unix/aganna/scripts/PCA/purcell5k.intervals

/home/unix/aganna/scripts/PCA/makeTPED.pl /broad/hptmp/aganna/outputPCA.vcf > /broad/hptmp/aganna/ouputPCA.tped

/home/unix/aganna/scripts/PCA/makeTFAM.pl /broad/hptmp/aganna/outputPCA.vcf > /broad/hptmp/aganna/outputPCA.tfam

/home/unix/aganna/scripts/PCA/combine_tped.pl /broad/hptmp/aganna/ouputPCA.tped > /broad/hptmp/aganna/MYOSEQ_1kg.tped

/home/unix/aganna/scripts/PCA/combine_tfam.pl /broad/hptmp/aganna/outputPCA.tfam > /broad/hptmp/aganna/MYOSEQ_1kg.tfam

plink --tfile /broad/hptmp/aganna/MYOSEQ_1kg --make-bed --out /broad/hptmp/aganna/MYOSEQ_1kg

/home/unix/aganna/scripts/PCA/prepare_input.pl /broad/hptmp/aganna/MYOSEQ_1kg

/home/unix/aganna/bin/smartpca -p /broad/hptmp/aganna/MYOSEQ_1kg.pca.par

/home/unix/aganna/scripts/PCA/prepare_plot.pl > ${outD}

rm /broad/hptmp/aganna/MYOSEQ*