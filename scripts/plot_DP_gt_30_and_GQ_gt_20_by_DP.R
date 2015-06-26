### PLOTS TI_TV by several VARIANTS characteristics  ###
### Argument 1: a file containing the % of variants with DP > 30 per sample
### Argument 2: a file containing the % of GQ > 20 by different DP levels per sample
### Argument 3: a folder where to output the plots

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
inD2 = as.character(args[2])
outD = as.character(args[3])

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/DPgt30"
#inD2="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/GQgt20byDP"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/plots/"

# Read data
d1DP30 <- fread(paste0(inD1,"_SNP_EX.txt"), header=F, stringsAsFactor=F)
d1GQ20byDP <- data.frame(fread(paste0(inD2,"_SNP_EX.txt"), header=F, stringsAsFactor=F))

d2DP30 <- fread(paste0(inD1,"_SNP_NEX.txt"), header=F, stringsAsFactor=F)
d2GQ20byDP <- data.frame(fread(paste0(inD2,"_SNP_NEX.txt"), header=F, stringsAsFactor=F))

d3DP30 <- fread(paste0(inD1,"_INDEL_EX.txt"), header=F, stringsAsFactor=F)
d3GQ20byDP <- data.frame(fread(paste0(inD2,"_INDEL_EX.txt"), header=F, stringsAsFactor=F))

d4DP30 <- fread(paste0(inD1,"_INDEL_NEX.txt"), header=F, stringsAsFactor=F)
d4GQ20byDP <- data.frame(fread(paste0(inD2,"_INDEL_NEX.txt"), header=F, stringsAsFactor=F))


d1GQ20byDP["V1"] <- NULL
d2GQ20byDP["V1"] <- NULL
d3GQ20byDP["V1"] <- NULL
d4GQ20byDP["V1"] <- NULL


newdata1DP30 <- data.frame(y=d1DP30$V2)
newdata1GQ20byDP <- data.frame(y=rowMeans(d1GQ20byDP,na.rm=T),x=seq(1,100))

newdata2DP30 <- data.frame(y=d2DP30$V2)
newdata2GQ20byDP <- data.frame(y=rowMeans(d2GQ20byDP,na.rm=T),x=seq(1,100))

newdata3DP30 <- data.frame(y=d3DP30$V2)
newdata3GQ20byDP <- data.frame(y=rowMeans(d3GQ20byDP,na.rm=T),x=seq(1,100))

newdata4DP30 <- data.frame(y=d4DP30$V2)
newdata4GQ20byDP <- data.frame(y=rowMeans(d4GQ20byDP,na.rm=T),x=seq(1,100))


png(paste0(outD,"DP_gt_30_by_sample.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y = y, x=1:nrow(newdata1DP30)), data=newdata1DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - SNPs, EXOMES") + ylim(0,1) + geom_hline(yintercept=mean(newdata1DP30$y),colour="yellow", linetype=1)

p2 <- ggplot(aes(y = y, x=1:nrow(newdata2DP30)), data=newdata2DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - SNPs, NON EXOMES") + ylim(0,1) + geom_hline(yintercept=mean(newdata2DP30$y),colour="yellow", linetype=1)

p3 <- ggplot(aes(y = y, x=1:nrow(newdata3DP30)), data=newdata3DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - INDELs, EXOMES") + ylim(0,1) + geom_hline(yintercept=mean(newdata3DP30$y),colour="yellow", linetype=1)

p4 <- ggplot(aes(y = y, x=1:nrow(newdata4DP30)), data=newdata4DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - INDELs, NON EXOMES") + ylim(0,1) + geom_hline(yintercept=mean(newdata4DP30$y),colour="yellow", linetype=1)

multiplot(p1, p2, p3, p4, cols = 2)

dev.off()


newdataall <- data.frame(y=c(rowMeans(d1GQ20byDP,na.rm=T),rowMeans(d2GQ20byDP,na.rm=T),rowMeans(d3GQ20byDP,na.rm=T),rowMeans(d4GQ20byDP,na.rm=T)),x=c(rep(seq(1,100),4)), group=factor(c(rep(1,100),rep(2,100),rep(3,100),rep(4,100))))

png(paste0(outD,"GQ_gt_20_by_DP.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y = y, x=x, group=group), data=newdataall) + geom_point(alpha=0.4) + stat_smooth(se = FALSE,lwd=2,span = 0.2, aes(colour=group)) + xlab("DP") + ylab("% of variants with GQ > 20") + xlim(c(0,100)) + ylim(c(0,1)) + ggtitle("GQ > 20 by DP") + scale_colour_discrete("Group",labels=c("SNPs-EXOMES","SNPs-NON EXOMES","INDELs EXOMES","INDELs-NON EXOMES"))

p1

dev.off()





