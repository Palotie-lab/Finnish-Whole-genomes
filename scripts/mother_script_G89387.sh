##########################
###### DEFINE PATHS ######
##########################

log="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/logs/G89387/"
plots="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/G89387/"
temp="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/"
measure="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/G89387/"
original="/humgen/atgu1/fs03/wip/aganna/fin_seq/original/seq/"
script="/humgen/atgu1/fs03/wip/aganna/fin_seq/scripts/seq/G89387/"


################################
#### RE-HEADER THE VCF FILE ####
################################

## Check the header ##
bcftools view -h ${original}G89387.vcf.gz > header.txt

#########################
#### BAM-LEVEL STATS ####
#########################

## Get GC content for each .bam file ##

for file in `awk -F $'\t' '{print $5}' /seq/dax/G89387/WGS/v12/G89387.calling_metadata.txt`; do
bsub -o ${log}gccontent/GC_`basename ${file} .bam`.log -R "rusage[mem=8]" -q week \
${script}calculate_GC_content.sh \
${file} \
${plots}gccontent \
${temp}gccontent;
done


## PLOT coverage_chimeras_insertsize_contamination ###
## and also save a file with per-sample coverage ###

Rscript ${script}plot_coverage_chimeras_insertsize_contamination.R \
${original}G89387.vcf.gz \
/seq/dax/G89387/WGS/v12/G89387.calling_metadata.txt \
${plots} ${measure}

## Plot GC content plots ##
## At the moment run it interactively ##
Rscript ${script}plot_GC_content.R \
${temp}gccontent/ \
${plots}G89387 \
${measure}


###################################################
#### POPULATION STRUCTURE: IBD, MDS, SEX CHECK ####
###################################################

## Run this script to convert the file in plink format ##
qsub -o ${log}to_plink_G89387.log -j y -b y -l m_mem_free=16g -q long -N TOPLINK -V \
/home/unix/aganna/scripts/vcf2plink.sh ${original}G89387.vcf.gz \
${original}G89387

## Run this script to overalp the .fam file with a correct .fam file ##
## Attention, notice that this script replace the original .fam file ##
#Rscript ${script}ann_fam.R \
#/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/G89387.fam \
#${original}pheno_gender.txt

## Calculate Kinship matrix, MDS and do a sex check ##
qsub -o ${log}kin_mds_sex.log -j y -b y -l m_mem_free=16g -q long -N KINMDS -V \
${script}calculate_kin_mds_sex_check.sh \
${original}G89387 \
${temp}G89387

## Plot stats  ##
Rscript ${script}plot_kin_mds.R \
${temp}G89387_QC2_pruned.kin0 \
${temp}G89387_QC2_prunedpc.ped \
${plots}



#######################################################
#### KEEP ONLY PASS AND EXCLUDE LCR AND BAD SAMPLES ###
#######################################################

## Keep PASS and remove low complexity regions ##
qsub -o ${log}keep_pass_G89387.log -j y -b y -l m_mem_free=16g -q long -N PASSNLC -V \
bcftools view -f .,PASS ${original}G89387.vcf.gz \
--force-samples \
-T ^/home/unix/aganna/LCR-hs37d5.bed \
-s ^MESTA_MEJ377,MESTA_MEJ308,MESTA_MEM125,MESTA_MET509 \
-Oz -o ${original}G89387_PASS_NLC.vcf.gz

tabix ${original}G89387_PASS_NLC.vcf.gz

## Plot general stats ##
bcftools stats ${original}G89387.vcf.gz > \
${measure}G89387.stats

bcftools stats ${original}G89387_PASS_NLC.vcf.gz > \
${measure}G89387_PASS_NLC.stats


##################################
#### OVERALAP PCs WITH EXAC ######
##################################

qsub -o ${log}G89387_pca_over_EXAC.log -j y -b y -l m_mem_free=8g -q long -N PCAEXAC -V \
${script}calculate_PCA_over_EXAC.sh \
${original}G89387_PASS_NLC.vcf.gz \
${temp}G89387_pca_over_EXAC.csv


######################################
#### PRODUCE DATASETS FOR ANALYSIS ###
######################################


## Will separate the original file in 4 parts, SNP/INDEL/EXOMES/WG and report statistics (other 3 additional files) ###
qsub -o ${log}G89387_SNP_EX.log -j y -b y -l m_mem_free=16g -q long -N SNP_EX -V \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G89387_PASS_NLC.vcf.gz \
${temp}G89387_PASS_NLC_SNP_EX.vcf "SNP" "YES" 

qsub -o ${log}G89387_SNP_WG.log -j y -b y -l m_mem_free=16g -q long -N SNP_WG -V \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G89387_PASS_NLC.vcf.gz \
${temp}G89387_PASS_NLC_SNP_WG.vcf "SNP" "NO" 

qsub -o ${log}G89387_INDEL_EX.log -j y -b y -l m_mem_free=16g -q long -N INDEL_EX -V \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G89387_PASS_NLC.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_EX.vcf "INDEL" "YES" 

qsub -o ${log}G89387_INDEL_WG.log -j y -b y -l m_mem_free=16g -q long -N INDEL_WG -V \
${script}calculate_EX_WG_SNP_INDEL.sh \
${original}G89387_PASS_NLC.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_WG.vcf "INDEL" "NO" 


## Save variants info in tabular format for plotting ##
## Obtain three additional file with genotype-specific values ##

qsub -o ${log}INFO_G89387_SNP_EX.log -j y -b y -l m_mem_free=16g -q long -N INFOSNPEX -V \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G89387_PASS_NLC_SNP_EX.vcf.gz \
${temp} G89387_PASS_NLC_SNP_EX

qsub -o ${log}INFO_G89387_SNP_WG.log -j y -b y -l m_mem_free=16g -q long -N INFOSNPWG -V \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G89387_PASS_NLC_SNP_WG.vcf.gz \
${temp} G89387_PASS_NLC_SNP_WG

qsub -o ${log}INFO_G89387_INDEL_EX.log -j y -b y -l m_mem_free=16g -q long -N INFOINDELEX -V \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G89387_PASS_NLC_INDEL_EX.vcf.gz \
${temp} G89387_PASS_NLC_INDEL_EX

qsub -o ${log}INFO_G89387_INDEL_WG.log -j y -b y -l m_mem_free=16g -q long -N INFOINDELWG -V \
${script}calculate_INFO_genospecific_AD_GT_DP_GQ.sh \
${temp}G89387_PASS_NLC_INDEL_WG.vcf.gz \
${temp} G89387_PASS_NLC_INDEL_WG


##### CHECK ####
## Calculate several metrics based on GQ and DP ##
bsub -o ${log}G89387_GQ_by_DP_SNP_EX.log  -R "rusage[mem=32]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G89387_PASS_NLC_SNP_EX.txt \
${temp}GQ_G89387_PASS_NLC_SNP_EX.txt \
${temp}samplenames_G89387_PASS_NLC_SNP_EX.txt \
${temp}

bsub -o ${log}G89387_GQ_by_DP_SNP_WG.log  -R "rusage[mem=90]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G89387_PASS_NLC_SNP_WG.txt \
${temp}GQ_G89387_PASS_NLC_SNP_WG.txt \
${temp}samplenames_G89387_PASS_NLC_SNP_WG.txt \
${temp}

bsub -o ${log}G89387_GQ_by_DP_INDEL_EX.log  -R "rusage[mem=16]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G89387_PASS_NLC_INDEL_EX.txt \
${temp}GQ_G89387_PASS_NLC_INDEL_EX.txt \
${temp}samplenames_G89387_PASS_NLC_INDEL_EX.txt \
${temp}

bsub -o ${log}G89387_GQ_by_DP_INDEL_WG.log  -R "rusage[mem=32]" -q priority \
python ${script}calculate_GQ_by_DP.py \
${temp}DP_G89387_PASS_NLC_INDEL_WG.txt \
${temp}GQ_G89387_PASS_NLC_INDEL_WG.txt \
${temp}samplenames_G89387_PASS_NLC_INDEL_WG.txt \
${temp}

####### CHECK ######


## Calculate AVERAGE GQ AND DP ##
qsub -o ${log}DP_GQ_DP_mean_G89387_SNP_EX.log -j y -b y -l m_mem_free=16g -q long -N DPGQMEANSNPEX -V \
${script}calculate_GQ_DP_mean.sh -f ${temp}G89387_PASS_NLC_SNP_EX.vcf.gz -o ${temp}DP_GQ_MEAN_G89387_PASS_NLC_SNP_EX.txt

qsub -o ${log}DP_GQ_DP_mean_G89387_SNP_WG.log -j y -b y -l m_mem_free=32g -q long -N DPGQMEANSNPWG -V \
${script}calculate_GQ_DP_mean.sh -f ${temp}G89387_PASS_NLC_SNP_WG.vcf.gz -o ${temp}DP_GQ_MEAN_G89387_PASS_NLC_SNP_WG.txt

qsub -o ${log}DP_GQ_DP_mean_G89387_INDEL_EX.log -j y -b y -l m_mem_free=16g -q long -N DPGQMEANINDELEX -V \
${script}calculate_GQ_DP_mean.sh -f ${temp}G89387_PASS_NLC_INDEL_EX.vcf.gz -o ${temp}DP_GQ_MEAN_G89387_PASS_NLC_INDEL_EX.txt

qsub -o ${log}DP_GQ_DP_mean_G89387_INDEL_WG.log -j y -b y -l m_mem_free=16g -q long -N DPGQMEANINDELWG -V \
${script}calculate_GQ_DP_mean.sh -f ${temp}G89387_PASS_NLC_INDEL_WG.vcf.gz -o ${temp}DP_GQ_MEAN_G89387_PASS_NLC_INDEL_WG.txt


## Calculate average indels lengths ##
qsub -o ${log}indel_length_G89387_INDEL_EX.log -j y -b y -l m_mem_free=8g -q long -N LENGTHEX -V \
${script}calculate_indel_length.sh -f ${temp}G89387_PASS_NLC_INDEL_EX.vcf.gz -o ${temp}indel_length_G89387_PASS_NLC_INDEL_EX.txt

qsub -o ${log}indel_length_G89387_INDEL_WG.log -j y -b y -l m_mem_free=8g -q long -N LENGTHWG -V \
${script}calculate_indel_length.sh -f ${temp}G89387_PASS_NLC_INDEL_WG.vcf.gz -o ${temp}indel_length_G89387_PASS_NLC_INDEL_WG.txt


## Convert each of the categories in  plink files ##
qsub -o ${log}G89387_plink_SNP_EX.log -j y -b y -l m_mem_free=32g -q long -N TOPLSNPEX -V \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G89387_PASS_NLC_SNP_EX.vcf.gz \
${temp}G89387_PASS_NLC_SNP_EX

qsub -o ${log}G89387_plink_SNP_WG.log -j y -b y -l m_mem_free=32g -q long -N TOPLSNPWG -V \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G89387_PASS_NLC_SNP_WG.vcf.gz \
${temp}G89387_PASS_NLC_SNP_WG

qsub -o ${log}G89387_plink_INDEL_EX.log -j y -b y -l m_mem_free=32g -q long -N TOPLINDELEX -V \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G89387_PASS_NLC_INDEL_EX.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_EX

qsub -o ${log}G89387_plink_INDEL_WG.log -j y -b y -l m_mem_free=32g -q long -N TOPLINDELWG -V \
/home/unix/aganna/scripts/vcf2plink.sh ${temp}G89387_PASS_NLC_INDEL_WG.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_WG


### Calculate stats for Heterozigosity, Hardy-weinberg equilibrium and Missing ###
qsub -o ${log}G89387_het_hardy_miss_SNP_EX.log -j y -b y -l m_mem_free=32g -q long -N SNPHHMEX -V \
${script}calculate_het_hardy_miss.sh \
${temp}G89387_PASS_NLC_SNP_EX

qsub -o ${log}G89387_het_hardy_miss_SNP_WG.log -j y -b y -l m_mem_free=32g -q long -N SNPHHMWG -V \
${script}calculate_het_hardy_miss.sh \
${temp}G89387_PASS_NLC_SNP_WG


qsub -o ${log}G89387_het_hardy_miss_INDEL_EX.log -j y -b y -l m_mem_free=32g -q long -N INDELHHMEX -V \
${script}calculate_het_hardy_miss.sh \
${temp}G89387_PASS_NLC_INDEL_EX

qsub -o ${log}G89387_het_hardy_miss_INDEL_WG.log -j y -b y -l m_mem_free=32g -q long -N INDELHHMWG -V \
${script}calculate_het_hardy_miss.sh \
${temp}G89387_PASS_NLC_INDEL_WG


## Run script on allele balance ##
qsub -o ${log}G89387_allele_balance_SNP_EX.log  -j y -b y -l m_mem_free=16g -q long -N ABSNPEX -V \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_SNP_EX.vcf.gz \
${temp}G89387_PASS_NLC_SNP_EX \
10000000

qsub -o ${log}G89387_allele_balance_SNP_WG.log  -j y -b y -l m_mem_free=16g -q long -N ABSNPWG -V \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_SNP_WG.vcf.gz \
${temp}G89387_PASS_NLC_SNP_WG \
10000000

qsub -o ${log}G89387_allele_balance_INDEL_EX.log  -j y -b y -l m_mem_free=16g -q long -N ABINDELEX -V \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_INDEL_EX.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_EX \
10000000

qsub -o ${log}G89387_allele_balance_INDEL_WG.log  -j y -b y -l m_mem_free=16g -q long -N ABINDELWG -V \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_INDEL_WG.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_WG \
10000000


###### TEMP ######
bsub -o ${log}G89387_allele_balance_SNP_EX.log  -R "rusage[mem=32]" -q priority \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_SNP_EX.vcf.gz \
${temp}G89387_PASS_NLC_SNP_EX \
10000000

bsub -o ${log}G89387_allele_balance_SNP_WG.log  -R "rusage[mem=32]" -q priority \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_SNP_WG.vcf.gz \
${temp}G89387_PASS_NLC_SNP_WG \
10000000

bsub -o ${log}G89387_allele_balance_INDEL_EX.log  -R "rusage[mem=32]" -q priority \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_INDEL_EX.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_EX \
10000000

bsub -o ${log}G89387_allele_balance_INDEL_WG.log  -R "rusage[mem=32]" -q priority \
python ${script}calculate_DP_80_20_z_score.py \
${temp}G89387_PASS_NLC_INDEL_WG.vcf.gz \
${temp}G89387_PASS_NLC_INDEL_WG \
10000000

####### TEMP ######


######################
#### PLOT RESULTS  ###
######################

# PCA over EXAC #
Rscript ${script}plot_pca_over_exac.R \
${temp}G89387_pca_over_EXAC.csv \
${plots}




## Plot stats from the GATK report for both the counteval and titveval evaluation modules ##
Rscript ${script}plot_gatk_stat.R \
${temp}COUNT_gatkreport_G89387_PASS_NLC \
${temp}TITV_gatkreport_G89387_PASS_NLC \
${measure}meanmedcoverage.csv \
${plots}

## Indels length ##
Rscript ${script}plot_indel_length.R \
${temp}indel_length_G89387_PASS_NLC \
${plots}


### Plot results from plink analysis ###
bsub -o ${log}plot_het_hardy_miss.log -R "rusage[mem=16]" -q priority \
Rscript ${script}plot_het_hardy_miss.R \
${temp}G77318RH \
${plots}



## Plot stats from bcftools report ##
Rscript ${script}plot_base_changes.R \
${temp}bcfreport_G77318RH \
${plots}


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
${temp}DPgt30_NLC \
${temp}GQgt20byDPgt_NLC \
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





#### KEEP ONLY GTs == 0 ####

bsub -o test.log -R "rusage[mem=8]" -q priority \
bcftools filter -i "FORMAT/GQ==0" --set-GTs "." /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_SNP_EX.vcf.gz -Oz -o G77318RH_SNP_EX_GQ0.vcf.gz

bsub -o test.log -R "rusage[mem=8]" -q priority \
bash test.txt

/home/unix/aganna/bcftools/bcftools filter  -i "FORMAT/GQ>0" --set-GTs "."  /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_SNP_EX.vcf.gz >  /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_EX_GQgt0.vcf
bgzip /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_EX_GQgt0.vcf
tabix /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_EX_GQgt0.vcf.gz




bcftools filter  -i "FORMAT/GQ>0" --set-GTs "."  /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_SNP_WG.vcf.gz >  /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_WG_GQgt0.vcf
bgzip /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_WG_GQgt0.vcf
tabix /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_WG_GQgt0.vcf.gz

bsub -o test3.log -R "rusage[mem=8]" -q priority \
${script}calculate_GQ_DP_mean.sh -f /humgen/atgu1/fs03/wip/aganna/G77318RH_SNP_WG_GQgt0.vcf.gz -o /humgen/atgu1/fs03/wip/aganna/DP_GQ_MEAN_G77318RH_SNP_WG_GQgt0.vcf.gz





qsub -o /humgen/atgu1/fs03/wip/aganna/test.log -j y -l m_mem_free=8g -q short -N test -V /humgen/atgu1/fs03/wip/aganna/test.txt


bsub -o test3.log -R "rusage[mem=8]" -q priority \
${script}calculate_GQ_DP_mean.sh -f G77318RH_SNP_EX_GQ0.vcf.gz -o DP_GQ_MEAN_G77318RH_SNP_EX_GQ0.vcf.gz

bsub -o test3.log -R "rusage[mem=8]" -q priority \
${script}calculate_GQ_DP_mean.sh -f G77318RH_SNP_EX_GQgt0.vcf.gz -o DP_GQ_MEAN_G77318RH_SNP_EX_GQgt0.vcf.gz


