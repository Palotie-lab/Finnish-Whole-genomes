##########################
###### DEFINE PATHS ######
##########################

log="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/logs/"
plots="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/"
temp="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/"
measure="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/"
original="/humgen/atgu1/fs03/wip/aganna/fin_seq/original/seq/"
script="/humgen/atgu1/fs03/wip/aganna/fin_seq/scripts/seq/"


################################
#### RE-HEADER THE VCF FILE ####
################################

## Add CPG island info ##
#java -jar /home/unix/aganna/snpEff/SnpSift.jar annotate -exists CPG /home/unix/aganna/cpg.vcf.gz ${original}G77318RH.vcf.gz | bgzip > /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RHCPG.vcf.gz  && tabix -f /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RHCPG.vcf.gz

## Save the header ##
bcftools view -h ${original}G77318.vcf.gz > ${temp}header.txt

## Replace all the 001 and 002 with _001 and _002 ##
sed -i '/#CHROM/c'"`cat ${temp}header.txt | grep "^#CHROM" | sed 's/\ 001/_001/g' | sed 's/\ 002/_002/g'`"'' ${temp}header.txt

## Reassign the new header to the .vcf file ##
## Latest version of bcftools is needed ##

/home/unix/aganna/bcftools/bcftools reheader -h ${temp}header.txt \
${original}G77318.vcf.gz \
-o ${original}G77318RH.vcf.gz
tabix ${original}G77318RH.vcf.gz

#########################
#### BAM-LEVEL STATS ####
#########################

## Get GC content for each .bam file ##

for file in `awk -F $'\t' '{print $5}' /seq/dax/G77318/WGS/v15/G77318.calling_metadata.txt`; do
bsub -o ${log}gccontent/GC_`basename ${file} .bam`.log -R "rusage[mem=8]" -q week \
${script}calculate_GC_content.sh \
${file} \
${plots}gccontent \
${temp}gccontent;
done


## PLOT coverage_chimeras_insertsize_contamination ###
## and also save a file with per-sample coverage ###

Rscript ${script}plot_coverage_chimeras_insertsize_contamination.R \
${original}G77318RH.vcf.gz \
/seq/dax/G77318/WGS/v15/G77318.calling_metadata.txt \
${plots} ${measure}

## Plot GC content plots ##
## At the moment run it interactively ##
Rscript ${script}plot_GC_content.R \
${temp}gccontent/ \
${plots} \
${measure}


###################################################
#### POPULATION STRUCTURE: IBD, MDS, SEX CHECK ####
###################################################

## Run this script to convert the file in plink format ##
bsub -o /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RH.log -R "rusage[mem=96]" -q week \
/home/unix/aganna/scripts/vcf2plink.sh ${original}G77318RH.vcf.gz \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RH

## Run this script to overalp the .fam file with a correct .fam file ##
## Attention, notice that this script replace the original .fam file ##
Rscript ${script}ann_fam.R \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RH.fam \
${original}pheno_gender.txt

## Calculate Kinship matrix, MDS and do a sex check ##
bsub -o ${temp}kin_mds_sex.log -R "rusage[mem=96]" -q week \
${script}calculate_kin_mds_sex_check.sh \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RH \
${temp}

## Plot stats  ##
Rscript ${script}plot_kin_mds.R \
${temp}G77318RH_QC1_pruned.kin0 \
${temp}G77318RH_QC1_prunedpc.ped \
${plots}



###################################
#### OBTAIN GENERAL STATISTICS ####
###################################

bcftools stats ${original}G77318RH.vcf.gz > \
${temp}G77318RH.stats


###############################
#### DECIDE VQSRT TRESHOLD ####
###############################

## Plot sample-wide stats ##
bsub -o ${log}INFO_G77318RH.log -R "rusage[mem=160]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${original}G77318RH.vcf.gz \
${temp} G77318RH

## Plot variants properties by VQSRTranche treshold ##
bsub -o ${log}INFO_G77318RH_plot.log -R "rusage[mem=72]" -q week \
Rscript ${script}plot_INFO_by_VQSRT_by_variant.R \
${temp}INFO_G77318RH.txt \
${plots} 597


### THESE ARE NEEDED TO SAVE QUANTITIES FOR FUTURE EXAMINATION ###

# Plink analysis #
bsub -o ${log}G77318RH_het_hardy_miss.log -R "rusage[mem=96]" -q week \
${script}calculate_het_hardy_miss.sh \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G77318RH


## Check concordance with chip data ##
bsub -o /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_EX.log  -R "rusage[mem=98]" -q week \
${script}calculate_concordance_WGS_chip.sh \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/FinnRisk \
${original}G77318RH.vcf.gz \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_G77318RH
  
## process GQ and DP ##  
bsub -o ${log}GQ_by_DP.log  -R "rusage[mem=230]" -q priority \
python ${script}calculate_GQ_by_DP2.py \
${temp}DP_G77318RH.txt \
${temp}GQ_G77318RH.txt \
${temp}samplenames_G77318RH_INDEL_WG.txt \
${temp}

## Allele balance ##
bsub -o ${log}allele_balance.log  -R "rusage[mem=64]" -q priority \
bash ${script}calculate_allele_balance.sh \
-f ${original}G77318RH.vcf.gz -o ${temp}AB_G77318RH.txt

## Run script for average DP and GQ ##
bsub -o ${log}DP_GQ_DP_mean.log -R "rusage[mem=64]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${original}G77318RH.vcf.gz -o${temp}DP_GQ_MEAN_G77318RH.txt


## Calculate allele balance and z scores from allele balance calculations ##
bsub -o ${log}allele_balance_python.log -R "rusage[mem=80]" -q priority \
python ${script}get_DP_80_20_z_score.py ${original}G77318RH.vcf.gz ${temp} 1000000

#bsub -o ${log}allele_balance_python2.log -R "rusage[mem=30]" -q priority \
#python ${script}get_DP_80_20_z_score.py ${temp}G77318RH_SNP_EX.vcf.gz ${temp}temp/ 100000

#python ${script}get_DP_80_20_z_score.py /humgen/atgu1/fs03/wip/aganna/fin_seq/original/playdata22.vcf.gz ${temp}temp/ 300



#######################
#### KEEP ONLY PASS ###
#######################

## Select samples to exclude ##
#d1 <- read.csv("/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/chimeras.csv", stringsAsFactor=F)
#d1$X.1[d1$chimeras > 0.05]

#d2 <- read.csv("/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/contamination.csv", stringsAsFactor=F)
#d1$X.1[d1$contamin > 0.05]

bcftools view -f .,PASS -S ${original}G77318RH.vcf.gz -Oz -o \
${original}G77318RH_PASS.vcf.gz



# Also do a sex check plot #

bcftools filter -r X ${original}G77318RH_PASS.vcf.gz -Oz -o chrX.vcf.gz
zcat chrX.vcf.gz | /home/unix/aganna/scripts/x_het.pl  > test.txt

zcat ${original}G77318RH_PASS.vcf.gz | /home/unix/aganna/scripts/y_cov.pl  > test2.txt


######################################
#### PRODUCE DATASETS FOR ANALYSIS ###
######################################


## Will separate the original file in 4 parts, SNP/INDEL/EXOMES/WG and report statistics (other 3 additional files) ###
bsub -o ${log}G77318RH_SNP_EX.log -R "rusage[mem=72]" -q week \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G77318RH_PASS.vcf.gz \
${temp}G77318RH_SNP_EX.vcf "SNP" "YES" 

bsub -o ${log}G77318RH_SNP_WG.log -R "rusage[mem=260]" -q week \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G77318RH_PASS.vcf.gz \
${temp}G77318RH_SNP_WG.vcf "SNP" "NO" 

bsub -o ${log}G77318RH_INDEL_EX.log -R "rusage[mem=72]" -q week \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G77318RH_PASS.vcf.gz \
${temp}G77318RH_INDEL_EX.vcf "INDEL" "YES" 

bsub -o ${log}G77318RH_INDEL_WG.log -R "rusage[mem=180]" -q week \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G77318RH_PASS.vcf.gz \
${temp}G77318RH_INDEL_WG.vcf "INDEL" "NO" 


## Save versions without low-complexity regions ##
bsub -o ${log}low_complexity_G77318RH_INDEL_EX.log -R "rusage[mem=15]" -q week \
bcftools view -T ^/home/unix/aganna/LCR-hs37d5.bed ${temp}G77318RH_INDEL_EX.vcf.gz -Oz -o ${temp}G77318RH_INDEL_EX_NLC.vcf.gz

bsub -o ${log}low_complexity_G77318RH_INDEL_WG.log -R "rusage[mem=15]" -q week \
bcftools view -T ^/home/unix/aganna/LCR-hs37d5.bed ${temp}G77318RH_INDEL_WG.vcf.gz -Oz -o ${temp}G77318RH_INDEL_WG_NLC.vcf.gz

bsub -o ${log}low_complexity_G77318RH_SNP_EX.log -R "rusage[mem=15]" -q week \
bcftools view -T ^/home/unix/aganna/LCR-hs37d5.bed ${temp}G77318RH_SNP_EX.vcf.gz -Oz -o ${temp}G77318RH_SNP_EX_NLC.vcf.gz

bsub -o ${log}low_complexity_G77318RH_SNP_WG.log -R "rusage[mem=15]" -q week \
bcftools view -T ^/home/unix/aganna/LCR-hs37d5.bed ${temp}G77318RH_SNP_WG.vcf.gz -Oz -o ${temp}G77318RH_SNP_WG_NLC.vcf.gz



## Save variants info in tabular format for plotting ##
## Obtain three additional file with genotype-specific values ##

bsub -o ${log}INFO_G77318RH_SNP_EX.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_SNP_EX.vcf.gz \
${temp} G77318RH_SNP_EX

bsub -o ${log}INFO_G77318RH_SNP_WG.log -R "rusage[mem=180]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_SNP_WG.vcf.gz \
${temp} G77318RH_SNP_WG

bsub -o ${log}INFO_G77318RH_INDEL_EX.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_INDEL_EX.vcf.gz \
${temp} G77318RH_INDEL_EX

bsub -o ${log}INFO_G77318RH_INDEL_WG.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_INDEL_WG.vcf.gz \
${temp} G77318RH_INDEL_WG


bsub -o ${log}INFO_G77318RH_INDEL_EX_NLC.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_INDEL_EX_NLC.vcf.gz \
${temp} G77318RH_INDEL_EX_NLC

bsub -o ${log}INFO_G77318RH_INDEL_WG_NLC.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_INDEL_WG_NLC.vcf.gz \
${temp} G77318RH_INDEL_WG_NLC



## Calculate several metrics based on GQ and DP ##
bsub -o ${log}GQ_by_DP_SNP_EX.log  -R "rusage[mem=98]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_SNP_EX.txt \
${temp}GQ_G77318RH_SNP_EX.txt \
${temp}samplenames_G77318RH_SNP_EX.txt \
${temp}

bsub -o ${log}GQ_by_DP_SNP_WG.log  -R "rusage[mem=200]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_SNP_WG.txt \
${temp}GQ_G77318RH_SNP_WG.txt \
${temp}samplenames_G77318RH_SNP_EX.txt \
${temp}

bsub -o ${log}GQ_by_DP_INDEL_EX.log  -R "rusage[mem=16]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_INDEL_EX.txt \
${temp}GQ_G77318RH_INDEL_EX.txt \
${temp}samplenames_G77318RH_INDEL_EX.txt \
${temp}

bsub -o ${log}GQ_by_DP_INDEL_WG.log  -R "rusage[mem=98]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_INDEL_WG.txt \
${temp}GQ_G77318RH_INDEL_WG.txt \
${temp}samplenames_G77318RH_INDEL_WG.txt \
${temp}


bsub -o ${log}GQ_by_DP_INDEL_EX_NLC.log  -R "rusage[mem=16]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_INDEL_EX_NLC.txt \
${temp}GQ_G77318RH_INDEL_EX_NLC.txt \
${temp}samplenames_G77318RH_INDEL_EX_NLC.txt \
${temp}

bsub -o ${log}GQ_by_DP_INDEL_WG_NLC.log  -R "rusage[mem=98]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_INDEL_WG_NLC.txt \
${temp}GQ_G77318RH_INDEL_WG_NLC.txt \
${temp}samplenames_G77318RH_INDEL_WG_NLC.txt \
${temp}





## Run script on allele balance ##
bash ${script}calculate_allele_balance.sh \
-f ${temp}G77318RH_SNP_EX.vcf.gz > \
${temp}AB_G77318RH_SNP_EX.txt

bash ${script}calculate_allele_balance.sh \
-f ${temp}G77318RH_SNP_WG.vcf.gz > \
${temp}AB_G77318RH_SNP_WG.txt


## Calculate AVERAGE GQ AND DP ##
bsub -o ${log}DP_GQ_DP_mean_SNP_EX.log -R "rusage[mem=16]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_SNP_EX.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_SNP_EX.txt

bsub -o ${log}DP_GQ_DP_mean_SNP_WG.log -R "rusage[mem=64]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_SNP_WG.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_SNP_WG.txt

bsub -o ${log}DP_GQ_DP_mean_INDEL_EX.log -R "rusage[mem=16]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_INDEL_EX.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_INDEL_EX.txt

bsub -o ${log}DP_GQ_DP_mean_INDEL_WG.log -R "rusage[mem=16]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_INDEL_WG.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_INDEL_WG.txt


bsub -o ${log}DP_GQ_DP_mean_SNP_EX_NLC.log -R "rusage[mem=16]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_SNP_EX_NLC.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_SNP_EX_NLC.txt

bsub -o ${log}DP_GQ_DP_mean_SNP_WG_NLC.log -R "rusage[mem=64]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_SNP_WG_NLC.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_SNP_WG_NLC.txt

bsub -o ${log}DP_GQ_DP_mean_INDEL_EX_NLC.log -R "rusage[mem=16]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_INDEL_EX_NLC.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_INDEL_EX_NLC.txt

bsub -o ${log}DP_GQ_DP_mean_INDEL_WG_NLC.log -R "rusage[mem=16]" -q priority \
${script}calculate_GQ_DP_mean.sh -f ${temp}G77318RH_INDEL_WG_NLC.vcf.gz -o ${temp}DP_GQ_MEAN_G77318RH_INDEL_WG_NLC.txt

## Calculate average indels lengths ##

bsub -o ${log}indel_length_INDEL_EX.log -R "rusage[mem=16]" -q priority \
${script}calculate_indel_length.sh -f ${temp}G77318RH_INDEL_EX.vcf.gz -o ${temp}indel_length_G77318RH_INDEL_EX.txt

bsub -o ${log}indel_length_INDEL_WG.log -R "rusage[mem=16]" -q priority \
${script}calculate_indel_length.sh -f ${temp}G77318RH_INDEL_WG.vcf.gz -o ${temp}indel_length_G77318RH_INDEL_WG.txt


bsub -o ${log}indel_length_INDEL_EX_NLC.log -R "rusage[mem=16]" -q priority \
${script}calculate_indel_length.sh -f ${temp}G77318RH_INDEL_EX_NLC.vcf.gz -o ${temp}indel_length_G77318RH_INDEL_EX_NLC.txt

bsub -o ${log}indel_length_INDEL_WG_NLC.log -R "rusage[mem=16]" -q priority \
${script}calculate_indel_length.sh -f ${temp}G77318RH_INDEL_WG_NLC.vcf.gz -o ${temp}indel_length_G77318RH_INDEL_WG_NLC.txt



## Convert each of the 4-categories files in  plink files ##
bsub -o ${log}G77318RH_plink_SNP_EX.log -R "rusage[mem=96]" -q week \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G77318RH_SNP_EX.vcf.gz \
${temp}G77318RH_SNP_EX

bsub -o ${log}G77318RH_plink_SNP_WG.log -R "rusage[mem=96]" -q week \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G77318RH_SNP_WG.vcf.gz \
${temp}G77318RH_SNP_WG

bsub -o ${log}G77318RH_plink_INDEL_EX.log -R "rusage[mem=96]" -q week \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G77318RH_INDEL_EX.vcf.gz \
${temp}G77318RH_INDEL_EX

bsub -o ${log}G77318RH_plink_INDEL_WG.log -R "rusage[mem=96]" -q week \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G77318RH_INDEL_WG.vcf.gz \
${temp}G77318RH_INDEL_WG


### Calculate stats for Heterozigosity, Hardy-weinberg equilibrium and Missing ###
bsub -o ${log}G77318RH_het_hardy_miss_SNP_EX.log -R "rusage[mem=96]" -q week \
${script}calculate_het_hardy_miss.sh \
${temp}G77318RH_SNP_EX

bsub -o ${log}G77318RH_het_hardy_miss_SNP_WG.log -R "rusage[mem=96]" -q week \
${script}calculate_het_hardy_miss.sh \
${temp}G77318RH_SNP_WG

bsub -o ${log}G77318RH_het_hardy_miss_INDEL_EX.log -R "rusage[mem=96]" -q week \
${script}calculate_het_hardy_miss.sh \
${temp}G77318RH_INDEL_EX

bsub -o ${log}G77318RH_het_hardy_miss_INDEL_WG.log -R "rusage[mem=96]" -q week \
${script}calculate_het_hardy_miss.sh \
${temp}G77318RH_INDEL_WG


######################
#### PLOT RESULTS  ###
######################


### Plot results from plink analysis ###
bsub -o ${log}plot_het_hardy_miss.log -R "rusage[mem=16]" -q priority \
Rscript ${script}plot_het_hardy_miss.R \
${temp}G77318RH \
${plots}


## Plot stats from the GATK report for both the counteval and titveval evaluation modules ##
Rscript ${script}plot_gatk_stat.R \
${temp}COUNT_gatkreport_G77318RH \
${temp}TITV_gatkreport_G77318RH \
${measure}meanmedcoverage.csv \
${plots}


## Plot stats from bcftools report ##
Rscript ${script}plot_base_changes.R \
${temp}bcfreport_G77318RH \
${plots}

## Indels length ##
Rscript ${script}plot_indel_length.R \
${temp}indel_length_G77318RH \
${plots}


## Plot some TI/TV by several variant-specific characteristics ##
Rscript ${script}plot_INFO_by_TI_TV_per_variant.R \
${temp}INFO_G77318RH \
${plots} 597

## Plot GQ and DP by variant ##
Rscript ${script}plot_GQ_DP_by_variant.R \
${temp}INFO_G77318RH \
${plots}
 

## Plot allele balance stats ##
Rscript ${script}plot_allele_balance_stat_by_variant.R \
${temp}AB_G77318RH \
${plots}
 
## Relationship between DP and GQ ##
Rscript ${script}plot_DP_GQ_relation.R \
${temp}DPgt30 \
${temp}GQgt20byDP \
${plots}
 
    
#################################################
##### CONCORDANCE ANALYSIS WITH PSYCH CHIP ######
#################################################

bsub -o /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_EX.log  -R "rusage[mem=98]" -q week \
${script}calculate_concordance_WGS_chip.sh \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/FinnRisk \
${temp}G77318RH_SNP_EX.vcf.gz \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_EX
    

bsub -o /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_WG.log  -R "rusage[mem=98]" -q week \
${script}calculate_concordance_WGS_chip.sh \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/FinnRisk \
${temp}G77318RH_SNP_WG.vcf.gz \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_WG


## Plot concordance ##
Rscript ${script}plot_concordance_WGS_chip.R \
/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance \
${plots}

 
########################################################
#### SCRIPT FOR CREATING BAF PLOTS FOR EVERY SAMPLE ####
########################################################

use .python-3.4.2 .r-3.2.1
unuse .gcc-4.9.0


/home/unix/giulio/socha/split.py --vcf ${original}G77318RH_PASS.vcf.gz --mask /psych/genetics_data/working/giulio/b37/beds/b37.clean.bed --out ${temp}bafplots/

for file in ${temp}bafplots/*.gz; 
    do /home/unix/giulio/socha/bafplot.py \
    --csv $file \
    --out ${plots}bafplots/`basename ${file} .gz`.png; \
    done
    
   
for i in `seq -f "%.0f" 1098411 1098458`;do
echo $i
bkill ${i};
done

for i in `seq 10 1 60`; do
bsub -o /broad/hptmp/aganna/${i}.log -R "rusage[mem=8]" -q hour \
${script}calculate_GQ20_DP.sh \
${temp}temp/G77318RH_INDEL_EX_ms.vcf.gz \
$i;
done

cat /broad/hptmp/aganna/GQ20DP_script_*_HET_temp  ${temp}temp/GQ20_DP_HET_EX.txt
cat /broad/hptmp/aganna/GQ20DP_script_*_ALL_temp  ${temp}temp/GQ20_DP_ALL_EX.txt

for i in `seq 10 1 60`; do
bsub -o /broad/hptmp/aganna/${i}.log -R "rusage[mem=8]" -q priority \
${script}calculate_GQ20_DP.sh \
${temp}temp/G77318RH_INDEL_WG_ms.vcf.gz \
$i;
done

cat /broad/hptmp/aganna/GQ20DP_script_*_HET_temp  ${temp}temp/GQ20_DP_HET_WG.txt
cat /broad/hptmp/aganna/GQ20DP_script_*_ALL_temp  ${temp}temp/GQ20_DP_ALL_WG.txt