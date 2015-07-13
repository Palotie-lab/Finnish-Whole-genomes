#!/usr/bin/env bash
###
# reports for all variants in vcf aan allelic balance deviation from 20/80 ratio
#  usage: get_allele_balance -f filename -n true -o outputfile
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

	bcftools query -f "%ID\t%DP\t%GQ_MEAN[\t%GT:%DP:%GQ]\n" $infile | awk 'BEGIN{ FS="\t"; OFS="\t"; print "ID","DP","GQ","N","N_GQLT20_DPGT30","N_NON_REF","DPTOT","DP_MEAN_NONREF","DEV_GQLT20_DPGT30"} { n_dev=0;  n_nonmissing=0; n_nonref=0; dpval=0; for(i=4; i<=NF;i++) { split($i,dat,":"); if(dat[1]!="./.") { n_nonmissing+=1; if(dat[2]>30 && dat[3]<20) n_dev+=1; if(dat[1]=="1/0"||dat[1]=="0/1"||dat[1]=="1/1") {n_nonref+=1; dpval+=dat[2]} } }  prop_dev=0; dp_dev=0; prop_dev=n_dev/n_nonmissing; if(n_nonref!=0) {dp_dev=dpval/n_nonref;} print $1,$2,$3,n_nonmissing,n_dev,n_nonref,dpval,dp_dev,prop_dev;  }' >$outfile

else
    
        $printer $infile | sed 's/^##FORMAT=<ID=AD,Number=\./##FORMAT=<ID=AD,Number=R/g' |  bcftools norm -Ou -m -any | bcftools norm -Ou -f $reference_seq | /home/unix/aganna/bcftools/bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' | bcftools query -f "%ID\t%DP\t%GQ_MEAN[\t%GT:%DP:%GQ]\n" | awk 'BEGIN{ FS="\t"; OFS="\t"; print "ID","DP","GQ","N_MISS","N_HET","N_HOMREF","N_HOMALT","DPTOT","DP_HET","DP_HOMREF","DP_HOMALT","GQTOT","GQ_HET","GQ_HOMREF","GQ_HOMALT"} {n_miss=0; n_het=0; n_homalt=0; n_homref=0; dp_het=0; gq_het=0; dp_homalt=0; gq_homalt=0; dp_homref=0; gq_homref=0; gqall=0; dpall=0; for(i=4; i<=NF;i++) { split($i,dat,":"); if (dat[1]=="./.") n_miss+=1; else { gqall+=dat[3]; dpall+=dat[2]; if(dat[1]=="1/0"||dat[1]=="0/1") {n_het+=1; dp_het+=dat[2]; gq_het+=dat[3];} else if (dat[1]=="1/1") {n_homalt+=1; dp_homalt+=dat[2]; gq_homalt+=dat[3];} else {n_homref+=1; dp_homref+=dat[2]; gq_homref+=dat[3];}}} print $1,$2,$3,n_miss,n_het,n_homref,n_homalt,dpall,dp_het,dp_homref,dp_homalt,gqall,gq_het,gq_homref,gq_homalt; }' >$outfile

fi
