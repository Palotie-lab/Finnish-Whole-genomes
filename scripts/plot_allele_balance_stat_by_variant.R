### ALLELE-BALANCE and % OF VARIANTS DEVIATION FROM 20/80 RATIO  ###
### Argument 1: a file containing the calculation needed (Alelle balance per variant, % of variants deviating the 80/20 ratio, GQ, AC)
### Argument 2: Directory location where the plots should be saved 

library(ggplot2)
library(data.table)
options(warn=-1)
library(reshape)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Read arguments #
args <- commandArgs(trailingOnly = TRUE)
inD1 = as.character(args[1])
outD = as.character(args[2])

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/AB_G77318RH"
#outD="/humgen/atgu1/fs03/wip/aganna/BDS_exomes/results/plots/"

inD1 = "/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/G89387_PASS_NLC"
inD2 = "/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/samplenames_G89387_PASS_NLC"
inD3 = "/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/DP_GQ_MEAN_G89387_PASS_NLC"


print("Reading data")
d1a <- fread(paste0(inD1,"_INDEL_EX_Zscores_per_sample.txt"), header=T, stringsAsFactor=F)
d1b <- fread(paste0(inD1,"_SNP_EX_Zscores_per_sample.txt"), header=T, stringsAsFactor=F)
d1c <- fread(paste0(inD1,"_INDEL_WG_Zscores_per_sample.txt"), header=T, stringsAsFactor=F)
d1d <- fread(paste0(inD1,"_SNP_WG_Zscores_per_sample.txt"), header=T, stringsAsFactor=F)

d2a <- fread(paste0(inD1,"_INDEL_EX_Albal_per_sample.txt"), header=T, stringsAsFactor=F)
d2b <- fread(paste0(inD1,"_SNP_EX_Albal_per_sample.txt"), header=T, stringsAsFactor=F)
d2c <- fread(paste0(inD1,"_INDEL_WG_Albal_per_sample.txt"), header=T, stringsAsFactor=F)
d2d <- fread(paste0(inD1,"_SNP_WG_Albal_per_sample.txt"), header=T, stringsAsFactor=F)

d3a <- read.csv(paste0(inD1,"_INDEL_EX_AB_stats_per_sample.txt"), header=T, stringsAsFactor=F)
d3b <- read.csv(paste0(inD1,"_SNP_EX_AB_stats_per_sample.txt"), header=T, stringsAsFactor=F)
d3c <- read.csv(paste0(inD1,"_INDEL_WG_AB_stats_per_sample.txt"), header=T, stringsAsFactor=F)
d3d <- read.csv(paste0(inD1,"_SNP_WG_AB_stats_per_sample.txt"), header=T, stringsAsFactor=F)


sampnamea <- read.table(paste0(inD2,"_INDEL_EX.txt"))
sampnameb <- read.table(paste0(inD2,"_SNP_EX.txt"))
sampnamec <- read.table(paste0(inD2,"_INDEL_WG.txt"))
sampnamed <- read.table(paste0(inD2,"_SNP_WG.txt"))


dgqdpa <- fread(paste0(inD3,"_INDEL_EX.txt"), header=T, stringsAsFactor=F)
dgqdpb <- fread(paste0(inD3,"_SNP_EX.txt"), header=T, stringsAsFactor=F)
dgqdpc <- fread(paste0(inD3,"_INDEL_WG.txt"), header=T, stringsAsFactor=F)
dgqdpd <- fread(paste0(inD3,"_SNP_WG.txt"), header=T, stringsAsFactor=F)


dgqdpa <- subset(dgqdpa, N_HET > 0)  
dgqdpb <- subset(dgqdpb, N_HET > 0)  
dgqdpc <- subset(dgqdpc, N_HET > 0)  
dgqdpd <- subset(dgqdpd, N_HET > 0)  



## Create unique ID ##
dgqdpa$POS <- sapply(strsplit(dgqdpa$ID,":"),"[[",2)
dgqdpa$CHR <- sapply(strsplit(dgqdpa$ID,":"),"[[",1)
dgqdpa$NEWID <- paste0(dgqdpa$CHR,":",dgqdpa$POS)

dgqdpb$POS <- sapply(strsplit(dgqdpb$ID,":"),"[[",2)
dgqdpb$CHR <- sapply(strsplit(dgqdpb$ID,":"),"[[",1)
dgqdpb$NEWID <- paste0(dgqdpb$CHR,":",dgqdpb$POS)

dgqdpc$POS <- sapply(strsplit(dgqdpc$ID,":"),"[[",2)
dgqdpc$CHR <- sapply(strsplit(dgqdpc$ID,":"),"[[",1)
dgqdpc$NEWID <- paste0(dgqdpc$CHR,":",dgqdpc$POS)

dgqdpd$POS <- sapply(strsplit(dgqdpd$ID,":"),"[[",2)
dgqdpd$CHR <- sapply(strsplit(dgqdpd$ID,":"),"[[",1)
dgqdpd$NEWID <- paste0(dgqdpd$CHR,":",dgqdpd$POS)


temp12a <- data.frame(NEWID=paste0(d2a[[1]],":",d2a[[2]]),AB=rowMeans(d2a[,3:ncol(d2a), with = FALSE], na.rm=T),BADAB=apply(d1a[,3:ncol(d1a), with = FALSE],1,function(x){sum(abs(x)>1.96,na.rm=T)/sum(!is.na(x),na.rm=T)}))
temp12b <- data.frame(NEWID=paste0(d2b[[1]],":",d2b[[2]]),AB=rowMeans(d2b[,3:ncol(d2b), with = FALSE], na.rm=T),BADAB=apply(d1b[,3:ncol(d1b), with = FALSE],1,function(x){sum(abs(x)>1.96,na.rm=T)/sum(!is.na(x),na.rm=T)}))

temp12c <- data.frame(NEWID=paste0(d2c[[1]],":",d2c[[2]]),AB=rowMeans(d2c[,3:ncol(d2c), with = FALSE], na.rm=T),BADAB=apply(d1c[,3:ncol(d1c), with = FALSE],1,function(x){sum(abs(x)>1.96,na.rm=T)/sum(!is.na(x),na.rm=T)}))

temp12d <- data.frame(NEWID=paste0(d2d[[1]],":",d2d[[2]]),AB=rowMeans(d2d[,3:ncol(d2d), with = FALSE], na.rm=T),BADAB=apply(d1d[,3:ncol(d1d), with = FALSE],1,function(x){sum(abs(x)>1.96,na.rm=T)/sum(!is.na(x),na.rm=T)}))



merge12a <- merge(dgqdpa,temp12a,by="NEWID")
merge12a <- merge12a[!duplicated(merge12a$NEWID) & merge12a$N_HET > 10,]
merge12b <- merge(dgqdpb,temp12b,by="NEWID")
merge12b <- merge12b[!duplicated(merge12b$NEWID) & merge12b$N_HET > 10,]

merge12c <- merge(dgqdpc,temp12c,by="NEWID")
merge12c <- merge12c[!duplicated(merge12c$NEWID) & merge12c$N_HET > 10,]

merge12d <- merge(dgqdpd,temp12d,by="NEWID")
merge12d <- merge12d[!duplicated(merge12d$NEWID) & merge12d$N_HET > 10,]


print("Now plotting")

png(paste0(outD,"allele_balance_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = AB), data=merge12b) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per HET variant, DP > 10, N HET > 10 - SNPs, EXOMES") + geom_vline(xintercept=mean(merge12b$AB, na.rm=T),colour="yellow", linetype=2) + geom_vline(xintercept=median(merge12b$AB, na.rm=T),colour="yellow", linetype=1)

p2 <- ggplot(aes(x = AB), data=merge12d) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per HET variant, DP > 10, N HET > 10 - SNPs, WG") + geom_vline(xintercept=mean(merge12d$AB, na.rm=T),colour="yellow", linetype=2) + geom_vline(xintercept=median(merge12d$AB, na.rm=T),colour="yellow", linetype=1)

p3 <- ggplot(aes(x = AB), data=merge12a) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per HET variant, DP > 10, N HET > 10 - INDELs, EXOMES") + geom_vline(xintercept=mean(merge12a$AB, na.rm=T),colour="yellow", linetype=2) + geom_vline(xintercept=median(merge12a$AB, na.rm=T),colour="yellow", linetype=1)

p4 <- ggplot(aes(x = AB), data=merge12c) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per HET variant, DP > 10, N HET > 10 - INDELs, WG") + geom_vline(xintercept=mean(merge12c$AB, na.rm=T),colour="yellow", linetype=2) + geom_vline(xintercept=median(merge12c$AB, na.rm=T),colour="yellow", linetype=1)


multiplot(p1, p2, p3, p4, cols = 2)

dev.off()



png(paste0(outD,"Bad_allele_balance_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = BADAB), data=merge12b) + geom_histogram() +ylab("Count") + xlab("% of samples with P(ALT|DP,0.46) < 0.05") + ggtitle(paste0("% of samples with bad AB per HET variant, DP > 10, N HET > 10 - SNPs, EXOMES \n N. of variants > 50% bad AB = ",sum(merge12b$BADAB > 0.5, na.rm=T)))

p2 <- ggplot(aes(x = BADAB), data=merge12d) + geom_histogram() +ylab("Count") + xlab("% of samples with P(ALT|DP,0.46) < 0.05") + ggtitle(paste0("% of samples with bad AB per HET variant, DP > 10, N HET > 10 - SNPs, WG \n N. of variants > 50% bad AB = ",sum(merge12d$BADAB > 0.5, na.rm=T)))

p3 <- ggplot(aes(x = BADAB), data=merge12a) + geom_histogram() +ylab("Count") + xlab("% of samples with P(ALT|DP,0.46) < 0.05") + ggtitle(paste0("% of samples with bad AB per HET variant, DP > 10, N HET > 10 - INDELs, EXOMES \n N. of variants > 50% bad ab = ",sum(merge12a$BADAB > 0.5, na.rm=T)))

p4 <- ggplot(aes(x = BADAB), data=merge12c) + geom_histogram() +ylab("Count") + xlab("% of samples with P(ALT|DP,0.46) < 0.05") + ggtitle(paste0("% of samples with bad AB per HET variant, DP > 10, N HET > 10 - INDELs, EXOMES \n N. of variants > 50% bad ab = ",sum(merge12c$BADAB > 0.5, na.rm=T)))

multiplot(p1, p2, p3, p4, cols = 2)

dev.off()


#head(merge12a[merge12a$BADAB==1,])

#bcftools query -r 10:104620163 /humgen/atgu1/fs03/wip/aganna/BDS_exomes/processed/seq/temp/C1836_PASS_NLC_INDEL.vcf.gz -f "%REF %ALT [%GT %SAMPLE %GQ %AD\n]"

#bcftools query -r 11:20648440 /humgen/atgu1/fs03/wip/aganna/BDS_exomes/processed/seq/temp/C1836_PASS_NLC_INDEL.vcf.gz -f "%REF %ALT [%GT %SAMPLE %GQ %AD\n]"

#bcftools query -r 9:78790216 /humgen/atgu1/fs03/wip/aganna/BDS_exomes/processed/seq/temp/C1836_PASS_NLC_INDEL.vcf.gz -f "%REF %ALT [%GT %SAMPLE %GQ %AD\n]"



print("Now plotting")

set.seed(123)
png(paste0(outD,"allele_balance_zscores_by_sample.png"), width=1000, height=1000, type="cairo")
par(mfrow=c(2,2))
plot (density(rnorm(100000)), col="red", lwd=5, xlim=c(-6,6), xlab="Z scores distribution for P(ALT|DP,0.46)", ylab="density", main="Z scores distribution P(ALT|DP,0.46) per sample, \n DP > 10, N HET > 10, SNPs-EXOMES")
for (i in 3:length(d1b))
{
	lines(density(d1b[[i]], na.rm=T),lwd=0.5)
	print(i)
}


plot (density(rnorm(100000)), col="red", lwd=5, xlim=c(-6,6), xlab="Z scores P(ALT|DP,0.46)", ylab="density", main="Z scores distribution P(ALT|DP,0.46) per sample, \n DP > 10, N HET > 10, SNPs-WG")
for (i in 3:length(d1d))
{
	lines(density(d1d[[i]], na.rm=T),lwd=0.5)
	print(i)
}
set.seed(123)
lines(density(rnorm(100000)), col="red", lwd=5)

set.seed(123)
lines(density(rnorm(100000)), col="red", lwd=5)


plot (density(rnorm(100000)), col="red", lwd=5, xlim=c(-6,6), xlab="Z scores P(ALT|DP,0.46)", ylab="density", main="Z scores distribution P(ALT|DP,0.46) per sample, \n DP > 10, N HET > 10, INDELs-EXOMES")
for (i in 3:length(d1a))
{
	lines(density(d1a[[i]], na.rm=T),lwd=0.5)
	print(i)
}
set.seed(123)
lines(density(rnorm(100000)), col="red", lwd=5)


plot (density(rnorm(100000)), col="red", lwd=5, xlim=c(-6,6), xlab="Z scores P(ALT|DP,0.46)", ylab="density", main="Z scores distribution P(ALT|DP,0.46) per sample, \n DP > 10, N HET > 10, INDELs-WG")
for (i in 3:length(d1c))
{
	lines(density(d1c[[i]], na.rm=T),lwd=0.5)
	print(i)
}
set.seed(123)
lines(density(rnorm(100000)), col="red", lwd=5)


dev.off()





# print("Structuring data")
#
# newdata1 <- data.frame(AB=as.numeric(d1$AB), ABP=as.numeric(d1$PROP_DEV_20_80)+0.001, AC=d1$AC, GQ=d1$GQ)
# newdata2 <- data.frame(AB=as.numeric(d2$AB), ABP=as.numeric(d2$PROP_DEV_20_80)+0.001, AC=d2$AC, GQ=d2$GQ)
#
#
# #nrow(newdata1[newdata1$AB < 0.2 | newdata1$AB > 0.8,])
# #nrow(newdata2[newdata2$AB < 0.2 | newdata2$AB > 0.8,])
#
#
#
# png(paste0(outD,"percent_bad_allele_balance_by_variant.png"), width=1200, height=800, type="cairo")
#
# p1 <- ggplot(aes(x = ABP), data=newdata1) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("% samples deviating from 20/80 allele balance per HET variant - SNPs, EXOMES") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))
#
# p1a <- ggplot(aes(x = ABP), data=newdata1[newdata1$AC>10,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("ONLY VARIANTS WITH AC > 10") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))
#
# p1b <- ggplot(aes(x = ABP), data=newdata1[newdata1$GQ>20,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("% samples deviating from 20/80 allele balance") + ggtitle("ONLY VARIANTS WITH GQ > 20") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))
#
#
# p2 <- ggplot(aes(x = ABP), data=newdata2) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("% samples deviating from 20/80 allele balance per HET variant - SNPs, WG") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))
#
# p2a <- ggplot(aes(x = ABP), data=newdata2[newdata2$AC>10,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("ONLY VARIANTS WITH AC > 10") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))
#
# p2b <- ggplot(aes(x = ABP), data=newdata2[newdata2$GQ>20,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("% samples deviating from 20/80 allele balance") + ggtitle("ONLY VARIANTS WITH GQ > 20") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))
#
#
# multiplot(p1, p1a,  p1b, p2, p2a,  p2b, cols = 2)
#
# dev.off()




