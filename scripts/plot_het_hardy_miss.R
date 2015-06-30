### VARIANT-SPECIFIC DP AND GQ PLOTS  ###
### Argument 1: the begining of file from a plink analysis which then gets read for the following extensions: hwe, het, imiss, lmiss
### Argument 2: Directory location where the plots should be saved 

library(ggplot2)
library(data.table)
options(warn=-1)

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
inD = as.character(args[1])
outD = as.character(args[2])

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/"
#NIND=597


print("Reading data")

d1a <- read.table(paste0(inD,"_SNP_EX.hwe"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d1b <- read.table(paste0(inD,"_SNP_EX.het"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d1c <- read.table(paste0(inD,"_SNP_EX.imiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d1d <- read.table(paste0(inD,"_SNP_EX.lmiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)

d2a <- read.table(paste0(inD,"_SNP_WG.hwe"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d2b <- read.table(paste0(inD,"_SNP_WG.het"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d2c <- read.table(paste0(inD,"_SNP_WG.imiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d2d <- read.table(paste0(inD,"_SNP_WG.lmiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)

d3a <- read.table(paste0(inD,"_INDEL_EX.hwe"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d3b <- read.table(paste0(inD,"_INDEL_EX.het"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d3c <- read.table(paste0(inD,"_INDEL_EX.imiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d3d <- read.table(paste0(inD,"_INDEL_EX.lmiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)

d4a <- read.table(paste0(inD,"_INDEL_WG.hwe"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d4b <- read.table(paste0(inD,"_INDEL_WG.het"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d4c <- read.table(paste0(inD,"_INDEL_WG.imiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)
d4d <- read.table(paste0(inD,"_INDEL_WG.lmiss"), header=T, stringsAsFactor=F, blank.lines.skip=F)



newdata1ac <- data.frame(HWE=as.numeric(d1a$P), MISSL=as.numeric(d1d$F_MISS))
newdata1bd <- data.frame(HET=as.numeric(d1b$F), MISSI=as.numeric(d1c$F_MISS), IDH=d1c$FID,IDM=d1c$FID)
newdata2ac <- data.frame(HWE=as.numeric(d2a$P), MISSL=as.numeric(d2d$F_MISS))
newdata2bd <- data.frame(HET=as.numeric(d2b$F), MISSI=as.numeric(d2c$F_MISS), IDH=d2c$FID,IDM=d2c$FID)
newdata3ac <- data.frame(HWE=as.numeric(d3a$P), MISSL=as.numeric(d3d$F_MISS))
newdata3bd <- data.frame(HET=as.numeric(d3b$F), MISSI=as.numeric(d3c$F_MISS), IDH=d3c$FID,IDM=d3c$FID)
newdata4ac <- data.frame(HWE=as.numeric(d4a$P), MISSL=as.numeric(d4d$F_MISS))
newdata4bd <- data.frame(HET=as.numeric(d4b$F), MISSI=as.numeric(d4c$F_MISS), IDH=d4c$FID,IDM=d4c$FID)

newdata1bd$IDH[!newdata1bd$IDH%in%newdata1bd$IDH[c(which.max(newdata1bd$HET),which.min(newdata1bd$HET))]] <- " "
newdata2bd$IDH[!newdata2bd$IDH%in%newdata2bd$IDH[c(which.max(newdata2bd$HET),which.min(newdata2bd$HET))]] <- " "
newdata3bd$IDH[!newdata3bd$IDH%in%newdata3bd$IDH[c(which.max(newdata3bd$HET),which.min(newdata3bd$HET))]] <- " "
newdata4bd$IDH[!newdata4bd$IDH%in%newdata4bd$IDH[c(which.max(newdata4bd$HET),which.min(newdata4bd$HET))]] <- " "

newdata1bd$IDM[newdata1bd$IDM!=newdata1bd$IDM[which.max(newdata1bd$MISSI)]] <- " "
newdata2bd$IDM[newdata2bd$IDM!=newdata2bd$IDM[which.max(newdata2bd$MISSI)]] <- " "
newdata3bd$IDM[newdata3bd$IDM!=newdata3bd$IDM[which.max(newdata3bd$MISSI)]] <- " "
newdata4bd$IDM[newdata4bd$IDM!=newdata4bd$IDM[which.max(newdata4bd$MISSI)]] <- " "



print("Now plotting")

png(paste0(outD,"HWE_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = HWE), data=newdata1ac) + 
    geom_bar(aes(y = (..count..)/sum(..count..))) + ylab("Proportion") + xlab("HWE P-value") + ggtitle(paste0("Hardy-Weinberg P-value per variant - SNPs, EXOMES \n N.of variants with P-value < 10^6=",sum(newdata1ac$HWE<0.000001), "(",round(sum(newdata1ac$HWE<0.000001)/length(newdata1ac$HWE),2),"%)")) + scale_x_log10() 

p2 <- ggplot(aes(x = HWE), data=newdata2ac) +  geom_bar(aes(y = (..count..)/sum(..count..))) +ylab("Proportion") + xlab("HWE P-value") + ggtitle(paste0("Hardy-Weinberg P-value per variant - SNPs, WG \n N.of variants with P-value < 10^6=",sum(newdata2ac$HWE<0.000001), "(",round(sum(newdata2ac$HWE<0.000001)/length(newdata2ac$HWE),2),"%)")) + scale_x_log10() 

p3 <- ggplot(aes(x = HWE), data=newdata3ac) +  geom_bar(aes(y = (..count..)/sum(..count..))) +ylab("Proportion") + xlab("HWE P-value") + ggtitle(paste0("Hardy-Weinberg P-value per variant - INDELs, EXOMES \n N.of variants with P-value < 10^6=",sum(newdata3ac$HWE<0.000001), "(",round(sum(newdata3ac$HWE<0.000001)/length(newdata3ac$HWE),2),"%)")) + scale_x_log10() 

p4 <- ggplot(aes(x = HWE), data=newdata4ac) +  geom_bar(aes(y = (..count..)/sum(..count..))) +ylab("Proportion") + xlab("HWE P-value") + ggtitle(paste0("Hardy-Weinberg P-value per variant - INDELs, WG \n N.of variants with P-value < 10^6=",sum(newdata4ac$HWE<0.000001), "(",round(sum(newdata4ac$HWE<0.000001)/length(newdata4ac$HWE),2),"%)")) + scale_x_log10() 

multiplot(p1, p2, p3, p4, cols = 2)

dev.off()



png(paste0(outD,"missing_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = MISSL), data=newdata1ac) + geom_histogram(origin = 0) + ylab("Count on log scale") + xlab("Missing proportion") + ggtitle(paste0("Proportion of missing per variant - SNPs, EXOMES \n N.of variants with missing > 5%=",sum(newdata1ac$MISSL>0.05), "(",round(sum(newdata1ac$MISSL>0.05)/length(newdata1ac$MISSL),5),"%)")) + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))


p2 <- ggplot(aes(x = MISSL), data=newdata2ac) + geom_histogram(origin = 0) + ylab("Count on log scale") + xlab("Missing proportion") + ggtitle(paste0("Proportion of missing per variant - SNPs, WG \n N.of variants with missing > 5%=",sum(newdata2ac$MISSL>0.05), "(",round(sum(newdata2ac$MISSL>0.05)/length(newdata2ac$MISSL),5),"%)")) + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))


p3 <- ggplot(aes(x = MISSL), data=newdata3ac) + geom_histogram(origin = 0) + ylab("Count on log scale") + xlab("Missing proportion") + ggtitle(paste0("Proportion of missing per variant - INDELs, EXOMES \n N.of variants with missing > 5%=",sum(newdata3ac$MISSL>0.05), "(",round(sum(newdata3ac$MISSL>0.05)/length(newdata3ac$MISSL),5),"%)")) + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))


p4 <- ggplot(aes(x = MISSL), data=newdata4ac) + geom_histogram(origin = 0) + ylab("Count on log scale") + xlab("Missing proportion") + ggtitle(paste0("Proportion of missing per variant - INDELs, WG \n N.of variants with missing > 5%=",sum(newdata4ac$MISSL>0.05), "(",round(sum(newdata4ac$MISSL>0.05)/length(newdata4ac$MISSL),5),"%)")) + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))


multiplot(p1, p2, p3, p4, cols = 2)

dev.off()




png(paste0(outD,"missing_by_sample.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y = MISSI, x=1:nrow(newdata1bd)), data=newdata1bd) + geom_point(size=4, col="red", alpha=0.8) + ylab("Missing proportion") + xlab("Sample") + ggtitle("Proportion of missing per sample - SNPs, EXOMEs") + geom_text(aes(label=IDM)) + geom_hline(yintercept=mean(newdata1bd$MISSI),colour="yellow", linetype=1)

p2 <- ggplot(aes(y = MISSI, x=1:nrow(newdata2bd)), data=newdata2bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("Missing proportion") + xlab("Sample") + ggtitle("Proportion of missing per sample - SNPs, WG") + geom_text(aes(label=IDM)) + geom_hline(yintercept=mean(newdata2bd$MISSI),colour="yellow", linetype=1)

p3 <- ggplot(aes(y = MISSI, x=1:nrow(newdata3bd)), data=newdata3bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("Missing proportion") + xlab("Sample") + ggtitle("Proportion of missing per sample - INDELs, EXOMEs") + geom_text(aes(label=IDM)) + geom_hline(yintercept=mean(newdata3bd$MISSI),colour="yellow", linetype=1)

p4 <- ggplot(aes(y = MISSI, x=1:nrow(newdata4bd)), data=newdata4bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("Missing proportion") + xlab("Sample") + ggtitle("Proportion of missing per sample - INDELs, WG") + geom_text(aes(label=IDM)) + geom_hline(yintercept=mean(newdata4bd$MISSI),colour="yellow", linetype=1)

multiplot(p1, p2, p3, p4, cols = 2)

dev.off()


png(paste0(outD,"F_by_sample.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y = HET, x=1:nrow(newdata1bd)), data=newdata1bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("F coef") + xlab("Sample") + ggtitle("F coefficient - SNPs, EXOMEs") + geom_text(aes(label=IDH)) + geom_hline(yintercept=mean(newdata1bd$HET),colour="yellow", linetype=1)

p2 <- ggplot(aes(y = HET, x=1:nrow(newdata1bd)), data=newdata2bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("F coef") + xlab("Sample") + ggtitle("F coefficient - SNPs, WG") + geom_text(aes(label=IDH)) + geom_hline(yintercept=mean(newdata2bd$HET),colour="yellow", linetype=1)

p3 <- ggplot(aes(y = HET, x=1:nrow(newdata1bd)), data=newdata3bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("F coef") + xlab("Sample") + ggtitle("F coefficient - IDNELs, EXOMEs") + geom_text(aes(label=IDH)) + geom_hline(yintercept=mean(newdata3bd$HET),colour="yellow", linetype=1)

p4 <- ggplot(aes(y = HET, x=1:nrow(newdata1bd)), data=newdata4bd) + geom_point(size=4, col="red", alpha=0.8) +ylab("F coef") + xlab("Sample") + ggtitle("F coefficient - INDELs, WG") + geom_text(aes(label=IDH)) + geom_hline(yintercept=mean(newdata4bd$HET),colour="yellow", linetype=1)

multiplot(p1, p2, p3, p4, cols = 2)

dev.off()







