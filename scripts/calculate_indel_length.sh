#!/usr/bin/env bash
###
# Report ID, REF and ALT allele and length of the indel
#  usage: calculate_indel_length -f filename -n true -o outputfile
#  optional argument -n true can be used if you have already split-multiallelics! This is essential for the script to 
#  work currently.
###


while getopts "f:n:o:" opt; do
	  declare "opt_$opt=${OPTARG:-0}"
done

infile=$opt_f
outfile=$opt_o

reference_seq="/humgen/1kg/reference/human_g1k_v37.fasta"

already_normalized=$opt_n

printer="cat"

if [ $infile=~"vcf.gz$" ]
then
	printer="zcat"
fi


if [ "$already_normalized" != "" ];
then

	bcftools query -f "%ID\t%REF\t%ALT\n" $infile | awk 'BEGIN{ FS="\t"; OFS="\t"; print "ID","REF","ALT","LENGTH"} { lenindel=0; lenindel=length($3)-length($2); print $1,$2,$3,lenindel; }' >$outfile

else
    
    $printer $infile | sed 's/^##FORMAT=<ID=AD,Number=\./##FORMAT=<ID=AD,Number=R/g' |  bcftools norm -Ou -m -any | bcftools norm -Ou -f $reference_seq | /home/unix/aganna/bcftools/bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' | bcftools query -f "%ID\t%REF\t%ALT\n" | awk 'BEGIN{ FS="\t"; OFS="\t"; print "ID","REF","ALT","LENGTH"} { lenindel=0; lenindel=length($3)-length($2); print $1,$2,$3,lenindel; }' >$outfile

fi
