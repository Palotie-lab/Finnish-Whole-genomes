### PLOT OVERLAP OF PCA WITH EXAC DATA ###
### Argument 1: the .csv file containing the PCs
### Argument 2: Directory location where the plots should be saved 

library(ggplot2)
options(warn=-1)


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

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

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/my_study_PCA.csv"
#outD="/humgen/atgu1/fs03/wip/aganna/BDS_exomes/results/plots/"

ExAC_PCA <- read.csv("/home/unix/aganna/scripts/PCA/ExAC_pca.csv",sep=",",header=T)
newstuff <- read.csv(inD,sep=",",header=T)

q2 <- ggplot(ExAC_PCA,aes(x=ExAC_PCA$PC1,y=ExAC_PCA$PC2)) + geom_point(alpha=0.10,size=0.85,color=ExAC_PCA$pop_color) + theme_bw()
q2 <- q2 + geom_point(data=newstuff,aes(x=newstuff$PC1,y=newstuff$PC2),alpha=0.80,size=1.2,color="black") 
q2 <- q2 + theme(legend.position="none")
q2 <- q2 + scale_x_continuous(name="PC1") + scale_y_continuous(name="PC2")

q3 <- ggplot(ExAC_PCA,aes(x=ExAC_PCA$PC2,y=ExAC_PCA$PC3)) + geom_point(alpha=0.05,size=0.85,color=ExAC_PCA$pop_color) + theme_bw()
q3 <- q3 + geom_point(data=newstuff,aes(x=newstuff$PC2,y=newstuff$PC3),alpha=0.80,size=1.2,color="black") 
q3 <- q3 + theme(legend.position="none")
q3 <- q3 + scale_x_continuous(name="PC2") + scale_y_continuous(name="PC3")

png(paste0(outD,"PCA_OVER_EXAC.png"),width=1200,height=600,type="cairo")
multiplot(q2,q3,cols=2)
dev.off()
