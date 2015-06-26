### PROGRAM TO PLOT COVERAGE, CHIMERIC %, INSERT SIZE, CONTAMINATION BY SAMPLE ###
### Argument 1: .vcf.gz file to analyze
### Argument 2: metadata file containing the location of .bam files
### Argument 3: dir where to output plots


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
origvcfD = as.character(args[1])
metadataD = as.character(args[2])
outD= as.character(args[3])
outDM = as.character(args[4])

#origvcfD="/humgen/atgu1/fs03/wip/aganna/fin_seq/original/seq/G77318RH.vcf.gz"
#metadataD="/seq/dax/G77318/WGS/v15/G77318.calling_metadata.txt"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/"
#outDM="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/"


# Extracting list of sampel ID from .vcf
a <- system(paste0('gzip -cd ',origvcfD, ' | grep  -m 1 "^#CHROM"'),intern = T)
b <- strsplit(paste0(a[1],a[2]),"\t")
c <- b[[1]]
c[length(c)] <- sub('/g', '', c[length(c)])
d <- c[10:length(c)]

# Read metadata #
metad <- read.table(metadataD, sep="\t", header=T, stringsAsFactor=F)

# Sample names have space, which is raplaced with _ #
metad$sample_name <- sub(" ","_",metad$sample_name)

# Keep only metadata also in the vcf ##
metads <- metad[metad$sample_name%in%d,]
metads$dir <- sapply(1:nrow(metads),function(x){sub(paste0(metads$sample_name[x],".bam"),"",metads$bam[x])})

# Function to extract data #
bametrics=function(metads,i)
{
	## Contamination ##
	adress <- paste0(metads$dir[i],metads$sample_name[i],".selfSM")
	ddt <- read.table(adress)
	contamin <- ddt$V8

	## Chimeras ##
	adress2 <- paste0(metads$dir[i],metads$sample_name[i],".alignment_summary_metrics")

	ddt2 <- read.table(adress2, comment.char = "#", sep="\t", header=T)
	chimeras <- ddt2$PCT_CHIMERAS[1]

	## Insert size ##
	adress3 <- paste0(metads$dir[i],metads$sample_name[i],".insert_size_metrics")

	ddt3 <- read.table(adress3, comment.char = "#", sep="\t", header=T, nrows = 2)
	insertmed <- ddt3$MEDIAN_INSERT_SIZE[1]

	insrtline <- as.numeric(strsplit(system(paste0('grep -nr -m 1 "^## HISTOGRAM" ',adress3),intern = T),":")[[1]][1])

	ddt4 <- read.table(adress3, comment.char = "#", sep="\t", header=T, skip = insrtline)

	insert <- ddt4[,2]

	## Coverage ##
	adress5 <- paste0(metads$dir[i],metads$sample_name[i],".wgs_metrics")

	ddt5 <- read.table(adress5, comment.char = "#", sep="\t", header=T,nrows = 1)

	meancoverage <- ddt5$MEAN_COVERAGE
	mediancoverage <- ddt5$MEDIAN_COVERAGE

	return(list(contamin,chimeras,meancoverage,mediancoverage,insert,insertmed))
}

## Get measures ##
bametricsR <- sapply(1:nrow(metads),bametrics,metads=metads)

## Process measures ##
contamin <- unlist(bametricsR[1,])	
chimeras <- unlist(bametricsR[2,])
meancoverage <- unlist(bametricsR[3,])
mediancoverage <- unlist(bametricsR[4,])
medianinsert <- unlist(bametricsR[6,])

write.csv(cbind(metads$sample_name,meancoverage,mediancoverage),file=paste0(outDM,"meanmedcoverage.csv"))
write.csv(cbind(metads$sample_name,chimeras),file=paste0(outDM,"chimeras.csv"))
write.csv(cbind(metads$sample_name,contamin),file=paste0(outDM,"contamination.csv"))
write.csv(cbind(metads$sample_name,medianinsert),file=paste0(outDM,"medianinsert.csv"))

## PLOT INSERT SIZE ###
inssizels <- NULL
for(i in 1:length(bametricsR[5,]))
{
	t <- cbind(metads$sample_name[i],bametricsR[5,i][[1]],seq(1:length(bametricsR[5,i][[1]])))
	inssizels <- rbind(inssizels,t)
}


newdata1 <- data.frame(group=inssizels[,1],y=as.numeric(inssizels[,2]),x=as.numeric(inssizels[,3]))
newdata1 <- newdata1[newdata1$x<600,]
meanins <- aggregate(newdata1$y,list(newdata1$x),median)


png(paste0(outD,"insert_size_distribution_by_sample.png"), width=800, height=800, type="cairo")
p1 <- ggplot(aes(y=y,x=x),data=newdata1) + geom_line(aes(group=group), alpha=0.5, colour="red") + xlab("Insert size") +ylab("Count")  + geom_line(aes(y=x,x=Group.1),data=meanins,colour="blue", lwd=2)
p1
dev.off()



### PLOT % CHIMERIC ###
newdata2 <- data.frame(y=chimeras,ID=metads$sample_name)
#newdata2 <- newdata2[order(newdata2$y),]

newdata2$ID[newdata2$ID!=newdata2$ID[which.max(newdata2$y)]] <- " "


png(paste0(outD,"chimeras_by_sample.png"), width=800, height=800, type="cairo")
p1 <- ggplot(aes(y = y, x=1:nrow(newdata2)), data=newdata2) + geom_point(size=10, alpha=0.5, colour="red") +xlab("Samples") + ylab("% Chimeric") + geom_text(aes(label=ID))
p1  
dev.off()


### PLOT COVERAGE ###
meancoverageN <- meancoverage[order(meancoverage)]
mediancoverageN <- mediancoverage[order(meancoverage)]

newdata3 <- data.frame(y=c(meancoverageN,mediancoverageN),x=c(seq(1:length(meancoverageN)),seq(1:length(mediancoverageN))),group=c(rep("Mean",length(meancoverageN)),rep("Median",length(mediancoverageN))),x=c(seq(1:length(meancoverageN))))

newdata3 <- newdata3[order(newdata3$y),]


png(paste0(outD,"coverage_by_sample.png"), width=1200, height=800, type="cairo")
p1 <- ggplot(aes(y = y, x=x), data=newdata3)  + geom_point(size=7, alpha=1) + geom_point(size=5, alpha=1,aes(colour=group)) +xlab("Samples") + ylab("Coverage") 
p1  
dev.off()


### PLOT % CONTAMINATION ###
newdata4 <- data.frame(y=contamin,ID=metads$sample_name)

newdata4$ID[newdata4$ID!=newdata4$ID[which.max(newdata4$y)]] <- " "

png(paste0(outD,"contamination_by_sample.png"), width=800, height=800, type="cairo")
p1 <- ggplot(aes(y = y, x=1:nrow(newdata4)), data=newdata4) + geom_point(size=10, alpha=0.5, colour="red") +xlab("Samples") + ylab("% Contamination") + geom_text(aes(label=ID))
p1  
dev.off()