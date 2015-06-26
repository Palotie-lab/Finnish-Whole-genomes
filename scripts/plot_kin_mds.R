### PROGRAM TO PLOT KINSHIP RELATIONSHIP ###
### Argument 1: a file generated from king (.kin0) including kinship relationship e.g. king -b file.bed --kinship
### Argument 2: a file generated from king (*pc.ped) including MDS e.g. king -b file.bed --mds
### Argument 3: folder where to save plots

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
inD1 = as.character(args[1])
inD2 = as.character(args[2])
outD = as.character(args[3])

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_QC1_pruned.kin0"
#inD2="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G77318RH_QC1_prunedpc.ped"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/plots/"

# Read data


d <- read.table(inD1, header=T, stringsAsFactor=F)

newdata <- data.frame(x=d$IBS0 , y=d$Kinship)
#newdata$y[newdata$y<0] <- 0
newdata$group <- ifelse(newdata$y > 0.354, "Twins/duplicates",
				 ifelse(newdata$y < 0.354 & newdata$y > 0.177, "1st degree",
				 ifelse(newdata$y < 0.177 & newdata$y > 0.0884, "2nd degree",
				 ifelse(newdata$y < 0.0884 & newdata$y > 0.0442, "3rd degree","Others"))))

png(paste0(outD,"IBD_vs_zero_IBS_sharing.png"), width=800, height=800, type="cairo")
ggplot(aes(x=x,y=y), data=newdata) + geom_point(aes(color=group)) + ggtitle("Kinship coefficient against the proportion of zero IBS-sharing") + ylab("Estimated Kinship Coefficient") + xlab("Proportion of Zero IBS")
dev.off()


## Plot MDS ##
d <- read.table(inD2, header=F, stringsAsFactor=F)

newdata1 <- data.frame(x=d$V7 , y=d$V8)
newdata2 <- data.frame(x=d$V7 , y=d$V9)
newdata3 <- data.frame(x=d$V8 , y=d$V9)
newdata4 <- data.frame(x=d$V9 , y=d$V10)

png(paste0(outD,"MDS.png"), width=800, height=800, type="cairo")
p1 <- ggplot(aes(x=x,y=y), data=newdata1) + geom_point(size=3, alpha=0.8) + ggtitle("MDS (distance 1 vs distance 2)") + xlab("Distance 1") + ylab("Distance 2")
p2 <- ggplot(aes(x=x,y=y), data=newdata2) + geom_point(size=3, alpha=0.8) + ggtitle("MDS (distance 1 vs distance 3)") + xlab("Distance 1") + ylab("Distance 3")
p3 <- ggplot(aes(x=x,y=y), data=newdata3) + geom_point(size=3, alpha=0.8) + ggtitle("MDS (distance 2 vs distance 3)") + xlab("Distance 2") + ylab("Distance 3")
p4 <- ggplot(aes(x=x,y=y), data=newdata4) + geom_point(size=3, alpha=0.8) + ggtitle("MDS (distance 3 vs distance 4)") + xlab("Distance 3") + ylab("Distance 4")

multiplot(p1,p2,p3,p4, cols=2)

dev.off()

