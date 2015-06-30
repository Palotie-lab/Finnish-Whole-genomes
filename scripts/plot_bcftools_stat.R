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

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/bcfreport_G77318RH"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/plots/"
# Read data


st1 <- system(paste0('grep "^ST" ', inD,'_SNP_EX.txt'), intern=TRUE)
st1m <- do.call(rbind, strsplit(st1,"\t"))

st2 <- system(paste0('grep "^ST" ', inD,'_SNP_WG.txt'), intern=TRUE)
st2m <- do.call(rbind, strsplit(st2,"\t"))

idd3 <- system(paste0('grep "^IDD" ', inD,'_INDEL_EX.txt'), intern=TRUE)
idd3m <- do.call(rbind, strsplit(idd3,"\t"))

idd4 <- system(paste0('grep "^IDD" ', inD,'_INDEL_WG.txt'), intern=TRUE)
idd4m <- do.call(rbind, strsplit(idd4,"\t"))




stdata1 <- data.frame(x=factor(st1m[,3]),y=as.numeric(st1m[,4]))
stdata2 <- data.frame(x=factor(st2m[,3]),y=as.numeric(st2m[,4]))

idddata3 <- data.frame(x=as.numeric(idd3m[,3]),y=as.numeric(idd3m[,4]))
idddata4 <- data.frame(x=as.numeric(idd4m[,3]),y=as.numeric(idd4m[,4]))

idddata3$lab=as.numeric(idddata3$x)
idddata3$lab[!idddata3$lab%in%c(3*seq(-60/3,60/3))] <- " "

idddata4$lab=as.numeric(idddata4$x)
idddata4$lab[!idddata4$lab%in%c(3*seq(-60/3,60/3))] <- " "


png(paste0(outD,"base_modifications.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y=y, x=x), data=stdata1) + geom_bar(stat="identity") +ylab("Count") + xlab("Type") + ggtitle("Base change - SNPs, EXOMES") +theme(axis.text=element_text(size=16,face="bold"))
p2 <- ggplot(aes(y=y, x=x), data=stdata2) + geom_bar(stat="identity") +ylab("Count") + xlab("Type") + ggtitle("Base change - SNPs, WG") +theme(axis.text=element_text(size=16,face="bold"))

multiplot(p1, p2, cols = 2)

dev.off()



png(paste0(outD,"length_indel.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y=y, x=x), data=idddata3) + geom_bar(stat="identity") +ylab("Count") + xlab("Length deletion - insertion") + ggtitle("Length indel - INDELs, EXOMES") +theme(axis.text=element_text(size=12,face="bold")) + geom_text(aes(label=lab), vjust=-1)
p2 <- ggplot(aes(y=y, x=x), data=idddata4) + geom_bar(stat="identity") +ylab("Count") + xlab("Length deletion - insertion") + ggtitle("Length indel - INDELs, WG") +theme(axis.text=element_text(size=12,face="bold")) + geom_text(aes(label=lab), vjust=-1)

multiplot(p1, p2, cols = 2)

dev.off()





