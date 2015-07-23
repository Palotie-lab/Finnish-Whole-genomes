### VARIANT-SPECIFIC DP AND GQ PLOTS  ###
### Argument 1: a file containing the requested fields for plotting (DP and GQ)
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
inD1 = as.character(args[1])
outD = as.character(args[2])
NIND = as.numeric(args[3])

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/DP_GQ_MEAN_G89387_PASS_NLC"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/G89387/"
#NIND=597

#inD1="/home/unix/aganna/"

#d1b <- fread(paste0(inD1,"DP_GQ_MEAN_G77318RH_SNP_EX_GQgt0.txt"), header=T, stringsAsFactor=F)

print("Reading data")
d1b <- fread(paste0(inD1,"_SNP_EX.txt"), header=T, stringsAsFactor=F)
d2b <- fread(paste0(inD1,"_SNP_WG.txt"), header=T, stringsAsFactor=F)
d3b <- fread(paste0(inD1,"_INDEL_EX.txt"), header=T, stringsAsFactor=F)
d4b <- fread(paste0(inD1,"_INDEL_WG.txt"), header=T, stringsAsFactor=F)



print("Structuring data")

newdata1 <- data.frame(DP=c(d1b$DP_HET[d1b$N_HET>0]/d1b$N_HET[d1b$N_HET>0],d1b$DP_HOMALT[d1b$N_HOMALT>0]/d1b$N_HOMALT[d1b$N_HOMALT>0],d1b$DP_HOMREF[d1b$N_HOMREF>0]/d1b$N_HOMREF[d1b$N_HOMREF>0]), group=c(rep("HET",sum(d1b$N_HET>0)),rep("HOM ALT",sum(d1b$N_HOMALT>0)),rep("HOM REF",sum(d1b$N_HOMREF>0))),GQ=c(d1b$GQ_HET[d1b$N_HET>0]/d1b$N_HET[d1b$N_HET>0],d1b$GQ_HOMALT[d1b$N_HOMALT>0]/d1b$N_HOMALT[d1b$N_HOMALT>0],d1b$GQ_HOMREF[d1b$N_HOMREF>0]/d1b$N_HOMREF[d1b$N_HOMREF>0]))

newdata2 <- data.frame(DP=c(d2b$DP_HET[d2b$N_HET>0]/d2b$N_HET[d2b$N_HET>0],d2b$DP_HOMALT[d2b$N_HOMALT>0]/d2b$N_HOMALT[d2b$N_HOMALT>0],d2b$DP_HOMREF[d2b$N_HOMREF>0]/d2b$N_HOMREF[d2b$N_HOMREF>0]), group=c(rep("HET",sum(d2b$N_HET>0)),rep("HOM ALT",sum(d2b$N_HOMALT>0)),rep("HOM REF",sum(d2b$N_HOMREF>0))),GQ=c(d2b$GQ_HET[d2b$N_HET>0]/d2b$N_HET[d2b$N_HET>0],d2b$GQ_HOMALT[d2b$N_HOMALT>0]/d2b$N_HOMALT[d2b$N_HOMALT>0],d2b$GQ_HOMREF[d2b$N_HOMREF>0]/d2b$N_HOMREF[d2b$N_HOMREF>0]))

newdata3 <- data.frame(DP=c(d3b$DP_HET[d3b$N_HET>0]/d3b$N_HET[d3b$N_HET>0],d3b$DP_HOMALT[d3b$N_HOMALT>0]/d3b$N_HOMALT[d3b$N_HOMALT>0],d3b$DP_HOMREF[d3b$N_HOMREF>0]/d3b$N_HOMREF[d3b$N_HOMREF>0]), group=c(rep("HET",sum(d3b$N_HET>0)),rep("HOM ALT",sum(d3b$N_HOMALT>0)),rep("HOM REF",sum(d3b$N_HOMREF>0))),GQ=c(d3b$GQ_HET[d3b$N_HET>0]/d3b$N_HET[d3b$N_HET>0],d3b$GQ_HOMALT[d3b$N_HOMALT>0]/d3b$N_HOMALT[d3b$N_HOMALT>0],d3b$GQ_HOMREF[d3b$N_HOMREF>0]/d3b$N_HOMREF[d3b$N_HOMREF>0]))

newdata4 <- data.frame(DP=c(d4b$DP_HET[d4b$N_HET>0]/d4b$N_HET[d4b$N_HET>0],d4b$DP_HOMALT[d4b$N_HOMALT>0]/d4b$N_HOMALT[d4b$N_HOMALT>0],d4b$DP_HOMREF[d4b$N_HOMREF>0]/d4b$N_HOMREF[d4b$N_HOMREF>0]), group=c(rep("HET",sum(d4b$N_HET>0)),rep("HOM ALT",sum(d4b$N_HOMALT>0)),rep("HOM REF",sum(d4b$N_HOMREF>0))),GQ=c(d4b$GQ_HET[d4b$N_HET>0]/d4b$N_HET[d4b$N_HET>0],d4b$GQ_HOMALT[d4b$N_HOMALT>0]/d4b$N_HOMALT[d4b$N_HOMALT>0],d4b$GQ_HOMREF[d4b$N_HOMREF>0]/d4b$N_HOMREF[d4b$N_HOMREF>0]))


#nrow(newdata1[newdata1$GQ<20,])
#nrow(newdata2[newdata2$GQ<20,])
#nrow(newdata3[newdata3$GQ<20,])
#nrow(newdata4[newdata4$GQ<20,])


print("Now plotting")


png(paste0(outD,"DP_by_variant.png"), width=1200, height=800, type="cairo")


p1 <- ggplot(aes(x = DP), data=newdata1) + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=3)+ ylab("Count") + xlab("Average Sequencing Depth (truncated to 60)") + ggtitle("Average Sequencing Depth per variant - SNPs, EXOMES") + xlim(c(0,60)) + geom_vline(xintercept=mean(newdata1$DP),colour="red", linetype=2)

p2 <- ggplot(aes(x = DP), data=newdata2)  + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=3)+ ylab("Count") + xlab("Average Sequencing Depth (truncated to 60)") + ggtitle("Average Sequencing Depth per variant - SNPs, WG") + xlim(c(0,60)) + geom_vline(xintercept=mean(newdata2$DP),colour="red", linetype=2)

p3 <- ggplot(aes(x = DP), data=newdata3)  + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=3)+ ylab("Count") + xlab("Average Sequencing Depth (truncated to 60)") + ggtitle("Average Sequencing Depth per variant - INDELs, EXOMES") + xlim(c(0,60)) + geom_vline(xintercept=mean(newdata3$DP),colour="red", linetype=2)

p4 <- ggplot(aes(x = DP), data=newdata4) + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=3)+ ylab("Count") + xlab("Average Sequencing Depth (truncated to 60)") + ggtitle("Average Sequencing Depth per variant - INDELs, WG") + xlim(c(0,60)) + geom_vline(xintercept=mean(newdata4$DP),colour="red", linetype=2)
multiplot(p1, p2, p3, p4, cols = 2)


dev.off()




png(paste0(outD,"GQ_MEAN_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = GQ), data=newdata1)  + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=6)+ ylab("Count") + xlab("Average Genotype quality") + ggtitle("Average Genotype quality per variant - SNPs, EXOMES")  + geom_vline(xintercept=mean(newdata1$GQ),colour="red", linetype=2)

p2 <- ggplot(aes(x = GQ), data=newdata2) + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=6)+ ylab("Count") + xlab("Average Genotype quality") + ggtitle("Average Genotype quality per variant - SNPs, WG")  + geom_vline(xintercept=mean(newdata2$GQ),colour="red", linetype=2)

p3 <- ggplot(aes(x = GQ), data=newdata3) + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=6)+ ylab("Count") + xlab("Average Genotype quality") + ggtitle("Average Genotype quality per variant - INDELs, EXOMES") + geom_vline(xintercept=mean(newdata3$GQ),colour="red", linetype=2)

p4 <- ggplot(aes(x = GQ), data=newdata4)  + geom_histogram(aes(fill=group),lwd=2,position="dodge", binwidth=6)+ ylab("Count") + xlab("Average Genotype quality") + ggtitle("Average Genotype quality per variant - INDELs, WG")  + geom_vline(xintercept=mean(newdata4$GQ),colour="red", linetype=2)

multiplot(p1, p2, p3, p4, cols = 2)

dev.off()


#newdata1$cat <- cut(newdata1$DP,seq(1,100,10))
#newdata2$cat <- cut(newdata2$DP,seq(1,100,10))
#newdata3$cat <- cut(newdata3$DP,seq(1,100,10))
#newdata4$cat <- cut(newdata4$DP,seq(1,100,10))


#aggregate(newdata1$GQ,list(newdata1$cat,newdata1$group),function(x){c(mean(x),length(x))})
#aggregate(newdata2$GQ,list(newdata2$cat,newdata2$group),function(x){c(mean(x),length(x))})
#aggregate(newdata3$GQ,list(newdata3$cat,newdata3$group),function(x){c(mean(x),length(x))})
#aggregate(newdata4$GQ,list(newdata4$cat,newdata4$group),function(x){c(mean(x),length(x))})

#d1b[d1b$DP_HOMREF/d1b$N_HOMREF > 50, ]


#bcftools query -r 2:905483 -f "[%GT %SAMPLE %AD %DP %GQ %PL\n]" G77318RH_SNP_EX_GQgt0.vcf

png(paste0(outD,"DP_by_GQ_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(data=newdata1,aes(y=GQ,x=DP)) + stat_binhex(colour="white",na.rm=TRUE, bin=50) + ylab("GQ - limited to 100") + xlab("DP - limited to 100") + xlim(c(0,100)) + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA, trans = "log",guide=FALSE)+ theme_bw() + ggtitle("DP by GQ - SNPs, EXOMES") + stat_smooth(aes(color=group),se = FALSE,lwd=3,span = 1) + ylim(c(0,100))

p2 <- ggplot(data=newdata2,aes(y=GQ,x=DP)) + stat_binhex(colour="white",na.rm=TRUE, bin=50) + ylab("GQ - limited to 100") + xlab("DP - limited to 100") + xlim(c(0,100)) + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA, trans = "log",guide=FALSE)+ theme_bw() + ggtitle("DP by GQ - SNPs, WG") + stat_smooth(aes(color=group),se = FALSE,lwd=3,span = 1) + ylim(c(0,100))

p3 <- ggplot(data=newdata3,aes(y=GQ,x=DP)) + stat_binhex(colour="white",na.rm=TRUE, bin=50) + ylab("GQ - limited to 100") + xlab("DP - limited to 100") + xlim(c(0,100)) + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA, trans = "log",guide=FALSE)+ theme_bw() + ggtitle("DP by GQ - INDELs, EXOMES") + stat_smooth(aes(color=group),se = FALSE,lwd=3,span = 1) + ylim(c(0,100))

p4 <- ggplot(data=newdata4,aes(y=GQ,x=DP)) + stat_binhex(colour="white",na.rm=TRUE, bin=50) + ylab("GQ - limited to 100") + xlab("DP - limited to 100") + xlim(c(0,100)) + scale_fill_gradientn(colours=c("green1","red"),name = "Frequency",na.value=NA, trans = "log",guide=FALSE)+ theme_bw() + ggtitle("DP by GQ - INDELs, WG") + stat_smooth(aes(color=group),se = FALSE,lwd=3,span = 1) + ylim(c(0,100))


multiplot(p1, p2, p3, p4, cols = 2)

dev.off()

