### PLOT STAT FROM GATK REPORT (TITVVARIANTSEVAL) ###
### Argument 1: a file generated from GATK VariantEval without the last part (which define each of the four dataset)
###				e.g. /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/TITV_gatkreport_study
### Argument 2: Directory location where the plots should be saved 

library(ggplot2)
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

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/TITV_gatkreport_G77318RH"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/plots/"

# Read data
d1 <- read.table(paste0(inD,"_SNP_EX.txt"), header=T, stringsAsFactor=F)
d2 <- read.table(paste0(inD,"_SNP_NEX.txt"), header=T, stringsAsFactor=F)

d1 <- d1[1:(nrow(d1)-1),]
d2 <- d2[1:(nrow(d2)-1),]


newdata1 <- data.frame(single=d1$tiTvRatio,IDtitv=as.character(d1$Sample))
newdata2 <- data.frame(single=d2$tiTvRatio,IDtitv=as.character(d2$Sample))

newdata1$IDtitv <- as.character(newdata1$IDtitv)
newdata2$IDtitv <- as.character(newdata2$IDtitv)

newdata1$IDtitv[!newdata1$IDtitv %in% newdata1$IDtitv[c(which.max(newdata1$single),which.min(newdata1$single))]] <- " "
newdata2$IDtitv[!newdata2$IDtitv %in% newdata2$IDtitv[c(which.max(newdata2$single),which.min(newdata2$single))]] <- " "


png(paste0(outD,"TITV_by_sample.png"), width=800, height=800, type="cairo")

p1 <- ggplot(aes(y = single, x=1:nrow(newdata1)), data=newdata1) + geom_point() +xlab("Samples") + ylab("Ti/Tv") + ggtitle("Ti/Tv - SNPs, EXOMES") + geom_text(aes(label=IDtitv))
p2 <- ggplot(aes(y = single, x=1:nrow(newdata2)), data=newdata2) + geom_point() +xlab("Samples") + ylab("Ti/Tv") + ggtitle("Ti/Tv - SNPs, NO EXOMES") + geom_text(aes(label=IDtitv))

multiplot(p1, p2, cols = 1)

dev.off()
