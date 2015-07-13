#!/usr/bin/env bash
### PROGRAM TO CALCULAT THE GX CONTENT FOR EACH .BAM FILE ###
### Argument 1: the .bam filename
### Argument 2: Folder where to output plots
### Argument 3: Folder where to output reports

inD="$1"
outD1="$2"
outD2="$3"

#inD=/seq/dax/G77318/WGS/v15/G77318.calling_metadata.txt
#outD1=${plots}gccontent
#outD2=${temp}gccontent

java -Xmx8g -jar /home/unix/aganna/picard-tools-1.135/picard.jar \
CollectGcBiasMetrics \
INPUT=${inD} \
CHART_OUTPUT=${outD1}/`basename ${inD} .bam`.pdf \
OUTPUT=${outD2}/`basename ${inD} .bam`.txt \
REFERENCE_SEQUENCE=/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta \
SUMMARY_OUTPUT=${outD2}/`basename ${inD} .bam`_summary.txt 


#-XX:ParallelGCThreads=4