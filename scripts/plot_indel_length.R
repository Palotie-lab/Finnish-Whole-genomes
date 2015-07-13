### PLOT STAT BCFTOOLS REPORT ###
### Argument 1: a file generated from bcftools stat without the last part (which define each of the four dataset)
###				e.g. /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/bcfreport_study
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

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/indel_length_G77318RH"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/"
# Read data


d1 <- read.table(paste0(inD,"_INDEL_EX.txt"), header=T, stringsAsFactor=F)
d2 <- read.table(paste0(inD,"_INDEL_WG.txt"), header=T, stringsAsFactor=F)

idddata1 <- data.frame(table(d1$LENGTH))
idddata2 <- data.frame(table(d2$LENGTH))

idddata1$lab=as.numeric(as.character(idddata1$Var1))
idddata1$lab[!idddata2$lab%in%c(3*seq(-60/3,60/3))] <- " "

idddata2$lab=as.numeric(as.character(idddata2$Var1))
idddata2$lab[!idddata2$lab%in%c(3*seq(-60/3,60/3))] <- " "


png(paste0(outD,"length_indel.png"), width=1200, height=800, type="cairo")

p3 <- ggplot(aes(y=Freq, x=as.numeric(as.character(Var1))), data=idddata1) + geom_bar(stat="identity") +ylab("Count") + xlab("Length deletion - insertion (truncated to 30)") + ggtitle("Length indel - INDELs, EXOMES - ZOOM ") +theme(axis.text=element_text(size=12,face="bold")) + geom_text(aes(label=lab), vjust=-1) + xlim(-30,30)
p4 <- ggplot(aes(y=Freq, x=as.numeric(as.character(Var1))), data=idddata2) + geom_bar(stat="identity") +ylab("Count") + xlab("Length deletion - insertion (truncated to 30)") + ggtitle("Length indel - INDELs, WG - ZOOM") +theme(axis.text=element_text(size=12,face="bold")) + geom_text(aes(label=lab), vjust=-1) + xlim(-30,30)


p1 <- ggplot(aes(y=Freq, x=as.numeric(as.character(Var1))), data=idddata1) + geom_bar(stat="identity") +ylab("Count") + xlab("Length deletion - insertion") + ggtitle("Length indel - INDELs, EXOMES") 
p2 <- ggplot(aes(y=Freq, x=as.numeric(as.character(Var1))), data=idddata2) + geom_bar(stat="identity") +ylab("Count") + xlab("Length deletion - insertion") + ggtitle("Length indel - INDELs, WG")  


multiplot(p1, p3, p2, p4, cols = 2)

dev.off()





