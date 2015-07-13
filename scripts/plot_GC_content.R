### PROGRAM TO PLOT COVERAGE, CHIMERIC %, INSERT SIZE, CONTAMINATION BY SAMPLE ###
### Argument 1: folder containing the results from the picard GC content analysis
### Argument 2: folder where to save the plots
### Argument 3: folder where to save the metrics


library(ggplot2)
options(warn=-1)
library(gridExtra)
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
outD1 = as.character(args[2])
outD2= as.character(args[3])


#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/gccontent/"
#outD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/"
#outD2="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/"

# Extracting list of sampel ID from .vcf
files <- list.files(inD)

ATDROP <- NULL
GCDROP <- NULL
WIN <- NULL
MEANDP <- NULL
BASQUAL <- NULL
FILE <- NULL
for (file in files)
{
	
	if (grepl("summary",file))
	{
		d <- read.table(paste0(inD,file),comment.char = "#", sep="\t", header=T)
		ATDROP <- c(ATDROP,d$AT_DROPOUT)
		GCDROP <- c(GCDROP,d$GC_DROPOUT)
	}
	else
	{
		d <- read.table(paste0(inD,file),comment.char = "#", sep="\t", header=T)
		WIN <- cbind(WIN,d$WINDOWS)
		MEANDP <- cbind(MEANDP,d$NORMALIZED_COVERAGE)
		BASQUAL <- cbind(BASQUAL,d$MEAN_BASE_QUALITY)
		FILE <- c(FILE,file)
	}
}


WINmean <- rowMeans(WIN)



newdata <- NULL
for(i in 1:ncol(MEANDP))
{
	newdata <- rbind(newdata,cbind(MEANDP[,i],1:101,FILE[i]))
}

newdata1 <- data.frame(newdata,stringsAsFactors=FALSE)
colnames(newdata1) <- c("y","x","group")
newdata1$x <- as.numeric(newdata1$x)
newdata1$y <- as.numeric(newdata1$y)

newdata2 <- data.frame(y=WINmean,x=1:101)

png(paste0(outD1,"GC_content_distribution_by_sample.png"), width=800, height=800, type="cairo")
p1 <- ggplot(aes(y=y,x=x),data=newdata1) + stat_smooth(aes(group=factor(group)), alpha=0.2, colour="red", se=F) + xlab("") +ylab("Fraction of normalized coverage (truncated at 2)")  + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm")) + xlim(c(1,90))

p2 <- ggplot(aes(y=y,x=x), data=newdata2) + geom_bar(stat="identity") + xlab("GC% of 100 base window") + ylab("Count") + ggtitle("") + theme( plot.margin=unit(c(-0.3,0.5,0.5,0), "cm")) + xlim(c(1,90))

grid.arrange(p1, p2, heights=c(2/3, 1/3), ncol=1)
dev.off()

newdata3 <- data.frame(y=c(ATDROP,GCDROP), x=c(1:(length(ATDROP)),1:(length(GCDROP))), group=c(rep("AT dropout",length(ATDROP)),rep("GC dropout",length(GCDROP))))
png(paste0(outD1,"AT_GC_dropout_by_sample.png"), width=800, height=800, type="cairo")
p1 <- ggplot(aes(y=y,x=x),data=newdata3) + geom_point(size=4, alpha=0.8, colour="red") + xlab("Samples") +ylab("%") + facet_wrap(~group, ncol=1, scales="free")
p1
dev.off()


#medsample <- apply(MEANDP,1,median)

#KS <- NULL
#for (k in 1:ncol(MEANDP))
#{
#	ks <- max(abs(medsample-d1DP30[[k]]))
#	KS <- c(KS,ks)
#}


toexp <- data.frame(at=ATDROP,gc=GCDROP,id=gsub(".txt","",FILE))
write.csv(toexp,paste0(outD2,"GC_measures.csv"), quote=F,row.names=F)