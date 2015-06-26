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

## PLOT coverage_chimeras_insertsize_contamination ###
## and also save a file with per-sample coverage ###

Rscript ${script}plot_coverage_chimeras_insertsize_contamination.R \
${original}G77318RH.vcf.gz \
/seq/dax/G77318/WGS/v15/G77318.calling_metadata.txt \
${plots} ${measure}


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

bsub -o ${log}INFO_G77318RH.log -R "rusage[mem=72]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${original}G77318RH.vcf.gz \
${temp} G77318RH

## Plot variants properties by VQSRTranche treshold ##
bsub -o ${log}INFO_G77318RH_plot.log -R "rusage[mem=72]" -q week \
Rscript ${script}plot_INFO_by_VQSRT_by_variant.R \
${temp}INFO_G77318RH.txt \
${plots} 597


#######################
#### KEEP ONLY PASS ###
#######################

bcftools view -f .,PASS ${original}G77318RH.vcf.gz -Oz -o \
${original}G77318RH_PASS.vcf.gz


######################################
#### PRODUCE DATASETS FOR ANALYSIS ###
######################################


## Will separate the original file in 4 parts, SNP/INDEL/EXOMES/WG and report statistics (other 3 additional files) ###
bsub -o ${log}G77318RH_SNP_EX.log -R "rusage[mem=72]" -q week \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G77318RH_PASS.vcf.gz \
${temp}G77318RH_SNP_EX.vcf "SNP" "YES" 

bsub -o ${log}G77318RH_SNP_WG.log -R "rusage[mem=180]" -q week \
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


## Save variants info in tabular format for plotting ##
## Obtain three additional file with genotype-specific values ##

bsub -o ${log}INFO_G77318RH_SNP_EX.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_SNP_EX.vcf.gz \
${temp} G77318RH_SNP_EX

bsub -o ${log}INFO_G77318RH_SNP_WG.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_SNP_WG.vcf.gz \
${temp} G77318RH_SNP_WG

bsub -o ${log}INFO_G77318RH_INDEL_EX.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_INDEL_EX.vcf.gz \
${temp} G77318RH_INDEL_EX

bsub -o ${log}INFO_G77318RH_INDEL_EX.log -R "rusage[mem=74]" -q week \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G77318RH_INDEL_WG.vcf.gz \
${temp} G77318RH_INDEL_WG


## Calculate several metrics based on GQ and DP ##
bsub -o ${log}GQ_by_DP_SNP_EX.log  -R "rusage[mem=98]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_SNP_EX.txt \
${temp}GQ_G77318RH_SNP_EX.txt \
${temp}

bsub -o ${log}GQ_by_DP_SNP_WG.log  -R "rusage[mem=200]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_SNP_WG.txt \
${temp}GQ_G77318RH_SNP_WG.txt \
${temp}

bsub -o ${log}GQ_by_DP_INDEL_EX.log  -R "rusage[mem=98]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_INDEL_EX.txt \
${temp}GQ_G77318RH_INDEL_EX.txt \
${temp}

bsub -o ${log}GQ_by_DP_INDEL_WG.log  -R "rusage[mem=98]" -q week \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G77318RH_INDEL_WG.txt \
${temp}GQ_G77318RH_INDEL_WG.txt \
${temp}


## Run script on allele balance ##
bash ${script}calculate_allele_balance.sh \
${temp}G77318RH_SNP_EX.vcf.gz -n false > \
${temp}AB_G77318RH_SNP_EX.txt

bash ${script}calculate_allele_balance.sh \
${temp}G77318RH_SNP_WG.vcf.gz -n false \
${temp}AB_G77318RH_SNP_WG.txt



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
Rscript ${script}plot_het_hardy_miss.R \
${temp}G77318RH \
${plots}


## Plot stats from the GATK report for both the counteval and titveval evaluation modules ##
Rscript ${script}plot_gatk_stat_COUNT.R \
${temp}COUNT_gatkreport_G77318RH \
${plots}

Rscript ${script}plot_gatk_stat_TI_TV.R \
${temp}TITV_gatkreport_G77318RH \
${plots}

## Plot stats from bcftools report ##
Rscript ${script}plot_bcftools_stat.R \
${temp}bcfreport_G77318RH \
${plots}


## Plot some TI/TV by several variant-specific characteristics ##
Rscript ${script}plot_INFO_by_TI_TV_per_variant.R \
${temp}INFO_G77318RH \
${plots} 597

## Plot GQ and DP by variant ##
Rscript ${script}plot_GQ_DP_by_variant.R \
${temp}INFO_G77318RH \
${plots} 597
 

## Plot allele balance stats ##
Rscript ${script}plot_allele_balance_stat_by_variant.R \
${temp}AB_G77318RH \
${temp}INFO_G77318RH \
${plots}
 

Rscript ${script}plot_DP_gt_30_and_GQ_gt_20_by_DP.R \
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