### ALLELE-BALANCE and % OF VARIANTS DEVIATION FROM 20/80 RATIO  ###
### Argument 1: a file containing the calculation needed (Alelle balance per variant, % of variants deviating the 80/20 ratio, GQ, AC)
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

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/AB_G77318RH"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/"


print("Reading data")
d1 <- fread(paste0(inD1,"_SNP_EX.txt"), header=T, stringsAsFactor=F)
d2 <- fread(paste0(inD1,"_SNP_WG.txt"), header=T, stringsAsFactor=F)


print("Structuring data")

newdata1 <- data.frame(AB=as.numeric(d1$AB), ABP=as.numeric(d1$PROP_DEV_20_80)+0.001, AC=d1$AC, GQ=d1$GQ)
newdata2 <- data.frame(AB=as.numeric(d2$AB), ABP=as.numeric(d2$PROP_DEV_20_80)+0.001, AC=d2$AC, GQ=d2$GQ)


print("Now plotting")

png(paste0(outD,"allele_balance_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = AB), data=newdata1) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per HET variant - SNPs, EXOMES") + geom_vline(xintercept=mean(newdata1$AB),colour="yellow", linetype=2) + geom_vline(xintercept=median(newdata1$AB),colour="yellow", linetype=1)

p2 <- ggplot(aes(x = AB), data=newdata2) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per HET variant - SNPs, WG") + geom_vline(xintercept=mean(newdata1$AB),colour="yellow", linetype=2) + geom_vline(xintercept=median(newdata1$AB),colour="yellow", linetype=1)

multiplot(p1, p2, cols = 2)

dev.off()


png(paste0(outD,"percent_bad_allele_balance_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = ABP), data=newdata1) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("% samples deviating from 20/80 allele balance per HET variant - SNPs, EXOMES") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))

p1a <- ggplot(aes(x = ABP), data=newdata1[newdata1$AC>10,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("ONLY VARIANTS WITH AC > 10") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))

p1b <- ggplot(aes(x = ABP), data=newdata1[newdata1$GQ>20,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("% samples deviating from 20/80 allele balance") + ggtitle("ONLY VARIANTS WITH GQ > 20") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))


p2 <- ggplot(aes(x = ABP), data=newdata2) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("% samples deviating from 20/80 allele balance per HET variant - SNPs, WG") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))

p2a <- ggplot(aes(x = ABP), data=newdata2[newdata2$AC>10,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("") + ggtitle("ONLY VARIANTS WITH AC > 10") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))

p2b <- ggplot(aes(x = ABP), data=newdata2[newdata2$GQ>20,]) + geom_histogram(origin = 0) +ylab("Count on log scale") + xlab("% samples deviating from 20/80 allele balance") + ggtitle("ONLY VARIANTS WITH GQ > 20") + scale_y_log10(breaks=c(10,100,1000,10000,100000,1000000,10000000),labels=c("10","100","1000","10000","100000","1000000","10000000"))


multiplot(p1, p1a,  p1b, p2, p2a,  p2b, cols = 2)

dev.off()




