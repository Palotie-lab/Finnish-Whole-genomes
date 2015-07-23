### PLOTS TI_TV by several VARIANTS characteristics  ###
### Argument 1: a file containing the % of variants with DP > 30 per sample
### Argument 2: a file containing the % of GQ > 20 by different DP levels per sample
### Argument 3: a folder where to output the plots

library(ggplot2)
library(data.table)
library(reshape)

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

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/DPgt30_NLC"
#inD2="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/GQgt20byDPgt_NLC"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/G89387/"

# Read data
d1DP30 <- fread(paste0(inD1,"_SNP_EX.txt"), header=F, stringsAsFactor=F)
d1GQ20byDPgt <- data.frame(fread(paste0(inD2,"_SNP_EX.txt"), header=T, stringsAsFactor=F))

d2DP30 <- fread(paste0(inD1,"_SNP_WG.txt"), header=F, stringsAsFactor=F)
d2GQ20byDPgt <- data.frame(fread(paste0(inD2,"_SNP_WG.txt"), header=T, stringsAsFactor=F))


d3DP30 <- fread(paste0(inD1,"_INDEL_EX.txt"), header=F, stringsAsFactor=F)
d3GQ20byDPgt <- data.frame(fread(paste0(inD2,"_INDEL_EX.txt"), header=T, stringsAsFactor=F))

d4DP30 <- fread(paste0(inD1,"_INDEL_WG.txt"), header=F, stringsAsFactor=F)
d4GQ20byDPgt <- data.frame(fread(paste0(inD2,"_INDEL_WG.txt"), header=T, stringsAsFactor=F))



d1GQ20byDPgt["V1"] <- NULL
d2GQ20byDPgt["V1"] <- NULL
d3GQ20byDPgt["V1"] <- NULL
d4GQ20byDPgt["V1"] <- NULL

newdata1DP30 <- data.frame(y=d1DP30$V2)
newdata2DP30 <- data.frame(y=d2DP30$V2)
newdata3DP30 <- data.frame(y=d3DP30$V2)
newdata4DP30 <- data.frame(y=d4DP30$V2)



png(paste0(outD,"DP_gt_30_by_sample.png"), width=1200, height=800, type="cairo")

p1 <- ggplot(aes(y = y, x=1:nrow(newdata1DP30)), data=newdata1DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - SNPs, EXOMES") + ylim(0,1) + geom_hline(yintercept=mean(newdata1DP30$y),colour="yellow", linetype=1)

p2 <- ggplot(aes(y = y, x=1:nrow(newdata2DP30)), data=newdata2DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - SNPs, WG") + ylim(0,1) + geom_hline(yintercept=mean(newdata2DP30$y),colour="yellow", linetype=1)

p3 <- ggplot(aes(y = y, x=1:nrow(newdata3DP30)), data=newdata3DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - INDELs, EXOMES") + ylim(0,1) + geom_hline(yintercept=mean(newdata3DP30$y),colour="yellow", linetype=1)

p4 <- ggplot(aes(y = y, x=1:nrow(newdata4DP30)), data=newdata4DP30) + geom_point(color="red", alpha=0.8, size=4) +xlab("Samples") + ylab("% of variants > 30x") + ggtitle("% of variants > 30x - INDELs, WG") + ylim(0,1) + geom_hline(yintercept=mean(newdata4DP30$y),colour="yellow", linetype=1)


multiplot(p1, p2, p3, p4, cols = 2)

dev.off()




d1GQ20byDPgt$mn <- rowMeans(d1GQ20byDPgt,na.rm=T)
d2GQ20byDPgt$mn <- rowMeans(d2GQ20byDPgt,na.rm=T)
d3GQ20byDPgt$mn <- rowMeans(d3GQ20byDPgt,na.rm=T)
d4GQ20byDPgt$mn <- rowMeans(d4GQ20byDPgt,na.rm=T)

meltex1 <- melt(d1GQ20byDPgt)
meltex2 <- melt(d2GQ20byDPgt)
meltex3 <- melt(d3GQ20byDPgt)
meltex4 <- melt(d4GQ20byDPgt)


meltex1$id <- paste0(meltex1$variable,"d1")
meltex2$id <- paste0(meltex2$variable,"d2")
meltex3$id <- paste0(meltex3$variable,"d3")
meltex4$id <- paste0(meltex4$variable,"d4")

meltex1$x <- rep(seq(5,75,5),length(d1GQ20byDPgt))
meltex2$x <- rep(seq(5,75,5),length(d2GQ20byDPgt))
meltex3$x <- rep(seq(5,75,5),length(d3GQ20byDPgt))
meltex4$x <- rep(seq(5,75,5),length(d4GQ20byDPgt))


meltex <- rbind(meltex1,meltex3)
meltex$min <- meltex$x-5

newdataall1 <- data.frame(y=meltex$value,x=meltex$x,min=meltex$min, id=factor(meltex$id), group=factor(c(rep(1,length(d1GQ20byDPgt)*nrow(d1GQ20byDPgt)),rep(3,length(d3GQ20byDPgt)*nrow(d3GQ20byDPgt)))))

newdataall1 <- newdataall1[!is.na(newdataall1$y),]
newdataall1$colorpoint <- factor(ifelse(newdataall1$id == "mnd3",2,
                      ifelse(newdataall1$id == "mnd1",1,NA)))
newdataall1$sizepoint <- ifelse(newdataall1$id %in% c("mnd1","mnd3"),8,0.001)
newdataall1$xnew <- ifelse(newdataall1$id %in% c("mnd1","mnd3"),newdataall1$x-2.5,NA)



meltex <- rbind(meltex2,meltex4)
meltex$min <- meltex$x-5

newdataall2 <- data.frame(y=meltex$value,x=meltex$x,min=meltex$min, id=factor(meltex$id), group=factor(c(rep(1,length(d2GQ20byDPgt)*nrow(d2GQ20byDPgt)),rep(3,length(d4GQ20byDPgt)*nrow(d4GQ20byDPgt)))))

newdataall2 <- newdataall2[!is.na(newdataall2$y),]
newdataall2$colorpoint <- factor(ifelse(newdataall2$id == "mnd4",2,
                      ifelse(newdataall2$id == "mnd2",1,NA)))
newdataall2$sizepoint <- ifelse(newdataall2$id %in% c("mnd2","mnd4"),8,0.001)
newdataall2$xnew <- ifelse(newdataall2$id %in% c("mnd2","mnd4"),newdataall2$x-2.5,NA)





png(paste0(outD,"GQ_gt_20_by_DP.png"), width=1200, height=800, type="cairo")

p1 <- ggplot( aes(x=y, ymin = min, ymax=x, group = id), data=newdataall1) + geom_linerange(aes(color=group), lwd=0.3, aplha=0.1) + geom_point(aes(y=xnew,size=sizepoint, fill=colorpoint), pch=21)   + xlab("% of variants with GQ > 20")  + xlim(c(0.5,1)) + ggtitle("% of GQ > 20 by DP bins - EXOMES") + theme(legend.text = element_text(size = 12)) + coord_flip() + scale_size_continuous(guide=FALSE)  + scale_y_continuous("DP Bins",breaks=seq(7.5,72.5,5),labels=c("5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74"),limits=c(5, 75)) + scale_fill_manual("",values=c("yellow","blue"),labels = c("Mean SNPs", "Mean INDELs")) + scale_colour_manual("Groups",values=c("yellow","blue"), labels = c("SNPs", "INDELs")) 

p2 <- ggplot( aes(x=y, ymin = min, ymax=x, group = id), data=newdataall2) + geom_linerange(aes(color=group), lwd=0.3, aplha=0.1) + geom_point(aes(y=xnew,size=sizepoint, fill=colorpoint), pch=21)   + xlab("% of variants with GQ > 20")  + xlim(c(0.5,1)) + ggtitle("% of GQ > 20 by DP bins - WG") + theme(legend.text = element_text(size = 12)) + coord_flip() + scale_size_continuous(guide=FALSE)  + scale_y_continuous("DP Bins",breaks=seq(7.5,72.5,5),labels=c("5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74"),limits=c(5, 75)) + scale_fill_manual("",values=c("yellow","blue"),labels = c("Mean SNPs", "Mean INDELs")) + scale_colour_manual("Groups",values=c("yellow","blue"), labels = c("SNPs", "INDELs")) 

multiplot(p1,p2,cols=1)

dev.off()






