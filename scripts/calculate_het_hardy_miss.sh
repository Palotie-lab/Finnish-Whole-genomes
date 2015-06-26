#!/bin/bash
### PROGRAM TO CALCULTE KINSHIP MATRIX, MDS AND SEX CHECK  ###
### Argument 1: Input file in plink (e.g. G77318RH). Notice the output use the input name and add the suffix for the specific analysis 

source /home/unix/aganna/.my.bashrc

inD="$1"
outD="$2"

### HETEROZYGOSITY ###
/home/unix/aganna/plink_linux_x86_64/plink -bfile ${inD} \
--het \
--out ${inD}

### H-W P-value ###
/home/unix/aganna/plink_linux_x86_64/plink -bfile ${inD} \
--hardy \
--out ${inD}

### % Missing data ###
/home/unix/aganna/plink_linux_x86_64/plink -bfile ${inD} \
--missing \
--out ${inD}

