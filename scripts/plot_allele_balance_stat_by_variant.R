### VARIANT-SPECIFIC DP AND GQ PLOTS  ###
### Argument 1: a file containing the allele balance per variant (just one column)
### Argument 2: a file containing the % samples deviating from 30/70 allele balance per each variant (just one column)
### Argument 3: Directory location where the plots should be saved 

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

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/AB"
#inD2="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/ABP"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/plots/"


print("Reading data")
d1a <- fread(paste0(inD1,"_SNP_EX.txt"), header=F, stringsAsFactor=F)
d1b <- fread(paste0(inD2,"_SNP_EX.txt"), header=F, stringsAsFactor=F)

d2a <- fread(paste0(inD1,"_SNP_NEX.txt"), header=F, stringsAsFactor=F)
d2b <- fread(paste0(inD2,"_SNP_NEX.txt"), header=F, stringsAsFactor=F)


print("Structuring data")

newdata1 <- data.frame(AB=as.numeric(d1a$V1), ABP=as.numeric(d1b$V1))
newdata2 <- data.frame(AB=as.numeric(d2a$V1), ABP=as.numeric(d2b$V1))


print("Now plotting")

png(paste0(outD,"allele_balance_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = AB), data=newdata1) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per variant - SNPs, EXOMES") + geom_vline(xintercept=mean(newdata1$AB),colour="yellow", linetype=2) + geom_vline(xintercept=median(newdata1$AB),colour="yellow", linetype=1)

p2 <- ggplot(aes(x = AB), data=newdata2) + geom_histogram() +ylab("Count") + xlab("Average allele balance") + ggtitle("Average allele balance per variant - SNPs, NON EXOMES") + geom_vline(xintercept=mean(newdata1$AB),colour="yellow", linetype=2) + geom_vline(xintercept=median(newdata1$AB),colour="yellow", linetype=1)

multiplot(p1, p2, cols = 2)

dev.off()


png(paste0(outD,"percent_bad_allele_balance_by_variant.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(x = ABP), data=newdata1) + geom_histogram() +ylab("Count") + xlab("% samples deviating from 30/70 allele balance") + ggtitle("% samples deviating from 30/70 allele balance per variant - SNPs, EXOMES") 

p2 <- ggplot(aes(x = ABP), data=newdata2) + geom_histogram() +ylab("Count") + xlab("% samples deviating from 30/70 allele balance") + ggtitle("% samples deviating from 30/70 allele balance per variant - SNPs, NON EXOMES") 

multiplot(p1, p2, cols = 2)

dev.off()


