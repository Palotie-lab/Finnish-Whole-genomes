### PLOT STAT FROM GATK REPORT  ###
### Argument 1: a file generated from GATK VariantEval without the last part (which define each of the four dataset)
###				e.g. /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/COUNT_gatkreport_study
### Argument 2: a file generated from GATK VariantEval without the last part (which define each of the four dataset)
###				e.g. /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/TITV_gatkreport_study
### Argument 3: location of the file including mean and median coverage
### Argument 4: Directory location where the plots should be saved 

library(ggplot2)
library(gridExtra)
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
coveD= as.character(args[3])
outD = as.character(args[4])

#inD1="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/COUNT_gatkreport_G89387_PASS_NLC"
#inD2="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/G89387/TITV_gatkreport_G89387_PASS_NLC"
#coveD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/measures/G89387/meanmedcoverage.csv"

#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/results/plots/G89387/"

# Read data (COUNT)
d1 <- read.table(paste0(inD1,"_SNP_EX.txt"), header=T, stringsAsFactor=F)
d2 <- read.table(paste0(inD1,"_SNP_WG.txt"), header=T, stringsAsFactor=F)
d3 <- read.table(paste0(inD1,"_INDEL_EX.txt"), header=T, stringsAsFactor=F)
d4 <- read.table(paste0(inD1,"_INDEL_WG.txt"), header=T, stringsAsFactor=F)

d1 <- d1[1:(nrow(d1)-1),]
d2 <- d2[1:(nrow(d2)-1),]
d3 <- d3[1:(nrow(d3)-1),]
d4 <- d4[1:(nrow(d4)-1),]

# Read data ti/tv ratio
# Read data
d1t <- read.table(paste0(inD2,"_SNP_EX.txt"), header=T, stringsAsFactor=F)
d2t <- read.table(paste0(inD2,"_SNP_WG.txt"), header=T, stringsAsFactor=F)

d1t <- d1t[1:(nrow(d1t)-1),]
d2t <- d2t[1:(nrow(d2t)-1),]

# Check that the samples are in the same order (they should be)
try(if (sum(d1$Sample==d1t$Sample) != nrow(d1)) stop("Not same samples in COUNT and TI/TV"))

newdata1 <- data.frame(single=d1$nSingletons , het_hom=d1$hetHomRatio, variant=d1$nVariantLoci, insdel=d1$insertionDeletionRatio, IDsing=as.character(d1$Sample), IDhet_hom=as.character(d1$Sample), IDvariant=as.character(d1$Sample), IDinsdel=as.character(d1$Sample),nInsertions=d1$nInsertions, nindel=d1$nInsertions+d1$nDeletions, nDeletions=d1$nDeletions, nHets=d1$nHets, nHom=d1$nHomRef + d1$nHomRef, titv=d1t$tiTvRatio,IDtitv=as.character(d1t$Sample))

newdata2 <- data.frame(single=d2$nSingletons , het_hom=d2$hetHomRatio, variant=d2$nVariantLoci, insdel=d2$insertionDeletionRatio, IDsing=as.character(d2$Sample), IDhet_hom=as.character(d2$Sample), IDvariant=as.character(d2$Sample), IDinsdel=as.character(d2$Sample),nInsertions=d2$nInsertions, nindel=d2$nInsertions+d2$nDeletions, nDeletions=d2$nDeletions, nHets=d2$nHets, nHom=d2$nHomRef + d2$nHomRef, titv=d2t$tiTvRatio,IDtitv=as.character(d2t$Sample))

newdata3 <- data.frame(single=d3$nSingletons , het_hom=d3$hetHomRatio, variant=d3$nVariantLoci, insdel=d3$insertionDeletionRatio, IDsing=as.character(d3$Sample), IDhet_hom=as.character(d3$Sample), IDvariant=as.character(d3$Sample), IDinsdel=as.character(d3$Sample),nInsertions=d3$nInsertions, nindel=d3$nInsertions+d3$nDeletions, nDeletions=d3$nDeletions, nHets=d3$nHets, nHom=d3$nHomRef + d3$nHomRef)

newdata4 <- data.frame(single=d4$nSingletons , het_hom=d4$hetHomRatio, variant=d4$nVariantLoci, insdel=d4$insertionDeletionRatio, IDsing=as.character(d4$Sample), IDhet_hom=as.character(d4$Sample), IDvariant=as.character(d4$Sample), IDinsdel=as.character(d4$Sample),nInsertions=d4$nInsertions, nindel=d4$nInsertions+d4$nDeletions, nDeletions=d4$nDeletions, nHets=d4$nHets, nHom=d4$nHomRef + d4$nHomRef)

newdata1$IDsing[!newdata1$IDsing %in% newdata1$IDsing[c(which.max(newdata1$single),which.min(newdata1$single))]] <- " "
newdata2$IDsing[!newdata2$IDsing %in% newdata2$IDsing[c(which.max(newdata2$single),which.min(newdata2$single))]] <- " "

newdata1$IDhet_hom[!newdata1$IDhet_hom %in% newdata1$IDhet_hom[c(which.max(newdata1$het_hom),which.min(newdata1$het_hom))]] <- " "
newdata2$IDhet_hom[!newdata2$IDhet_hom %in% newdata2$IDhet_hom[c(which.max(newdata2$het_hom),which.min(newdata2$het_hom))]] <- " "
newdata3$IDhet_hom[!newdata3$IDhet_hom %in% newdata3$IDhet_hom[c(which.max(newdata3$het_hom),which.min(newdata3$het_hom))]] <- " "
newdata4$IDhet_hom[!newdata4$IDhet_hom %in% newdata4$IDhet_hom[c(which.max(newdata4$het_hom),which.min(newdata4$het_hom))]] <- " "

newdata1$IDvariant[!newdata1$IDvariant %in% newdata1$IDvariant[c(which.max(newdata1$variant),which.min(newdata1$variant))]] <- " "
newdata2$IDvariant[!newdata2$IDvariant %in% newdata2$IDvariant[c(which.max(newdata2$variant),which.min(newdata2$variant))]] <- " "
newdata3$IDvariant[!newdata3$IDvariant %in% newdata3$IDvariant[c(which.max(newdata3$variant),which.min(newdata3$variant))]] <- " "
newdata4$IDvariant[!newdata4$IDvariant %in% newdata4$IDvariant[c(which.max(newdata4$variant),which.min(newdata4$variant))]] <- " "

newdata3$IDinsdel[!newdata3$IDinsdel %in% newdata3$IDinsdel[c(which.max(newdata3$insdel),which.min(newdata3$insdel))]] <- " "
newdata4$IDinsdel[!newdata4$IDinsdel %in% newdata4$IDinsdel[c(which.max(newdata4$insdel),which.min(newdata4$insdel))]] <- " "

newdata1$IDtitv <- as.character(newdata1$IDtitv)
newdata2$IDtitv <- as.character(newdata2$IDtitv)

newdata1$IDtitv[!newdata1$IDtitv %in% newdata1$IDtitv[c(which.max(newdata1$titv),which.min(newdata1$titv))]] <- " "
newdata2$IDtitv[!newdata2$IDtitv %in% newdata2$IDtitv[c(which.max(newdata2$titv),which.min(newdata2$titv))]] <- " "



DPd <- read.csv(coveD)
colnames(DPd) <- c("IND","ID","mean","median")

newdata1DP <- merge(newdata1,DPd,by.x="IDinsdel",by.y="ID")
newdata1DP <- newdata1DP[order(newdata1DP$mean), ]

newdata2DP <- merge(newdata2,DPd,by.x="IDinsdel",by.y="ID")
newdata2DP <- newdata2DP[order(newdata2DP$mean), ]


png(paste0(outD,"singeltons_by_sample.png"), width=800, height=800, type="cairo")

p1 <- ggplot(aes(y = single, x=1:nrow(newdata1DP)), data=newdata1DP) + geom_point(size=4,alpha=0.8,aes(color=mean)) +xlab("Samples ordered by DP") + ylab("N. Singlentons") + ggtitle("N. Singletons - SNPs, EXOMES") + geom_text(aes(label=IDsing)) + scale_colour_gradient("Mean coverage",low = "blue", high="red")
p2 <- ggplot(aes(y = single, x=1:nrow(newdata2DP)), data=newdata2DP) + geom_point(size=4,alpha=0.8,aes(color=mean)) + scale_colour_gradient("Mean coverage",low = "blue", high="red") +xlab("Samples ordered by DP") + ylab("N. Singlentons") + ggtitle("N. Singletons - SNPs, WG") + geom_text(aes(label=IDsing))

multiplot(p1, p2, cols = 1)

dev.off()


newdata1o <- newdata1[order(newdata1$variant),]
newdata2o <- newdata1[order(newdata2$variant),]
newdata3o <- newdata1[order(newdata3$variant),]
newdata4o <- newdata1[order(newdata4$variant),]


png(paste0(outD,"hetHomRatio_by_sample.png"), width=1600, height=800, type="cairo")

p1 <- ggplot(aes(y = het_hom, x=1:nrow(newdata1o)), data=newdata1o) + geom_point(size=4,alpha=0.8) +xlab("Samples") + ylab("Het/Hom") + ggtitle("Het/Hom - SNPs, EXOMES") + geom_text(aes(label=IDhet_hom)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))
p1a <- ggplot(aes(y = nHets, x=1:nrow(newdata1o)), data=newdata1o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. heterozygous")  + xlim(quantile(d2$VQSLOD,0.1),quantile(d2$VQSLOD,0.9)) + theme( plot.margin=unit(c(0.5,0.5,-1,0), "cm"))
p1b <- ggplot(aes(y = nHom, x=1:nrow(newdata1o)), data=newdata1o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. homozygous")  + xlim(quantile(d2$VQSLOD,0.1),quantile(d2$VQSLOD,0.9)) + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm")) + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))

p2 <- ggplot(aes(y = het_hom, x=1:nrow(newdata2o)), data=newdata2o) + geom_point(size=4,alpha=0.8) +xlab("Samples") + ylab("Het/Hom") + ggtitle("Het/Hom - SNPs, WG") + geom_text(aes(label=IDhet_hom)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))
p2a <- ggplot(aes(y = nHets, x=1:nrow(newdata2o)), data=newdata2o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. heterozygous") + theme( plot.margin=unit(c(0.5,0.5,-1,0), "cm"))
p2b <- ggplot(aes(y = nHom, x=1:nrow(newdata2o)), data=newdata2o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. homozygous") + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))


p3 <- ggplot(aes(y = het_hom, x=1:nrow(newdata3o)), data=newdata3o) + geom_point(size=4,alpha=0.8) +xlab("Samples") + ylab("Het/Hom") + ggtitle("Het/Hom - INDELSs, EXOMES") + geom_text(aes(label=IDhet_hom)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))
p3a <- ggplot(aes(y = nHets, x=1:nrow(newdata3o)), data=newdata3o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. heterozygous") + theme( plot.margin=unit(c(0.5,0.5,-1,0), "cm"))
p3b <- ggplot(aes(y = nHom, x=1:nrow(newdata3o)), data=newdata3o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. homozygous") + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))


p4 <- ggplot(aes(y = het_hom, x=1:nrow(newdata4o)), data=newdata4o) + geom_point(size=4,alpha=0.8) +xlab("Samples") + ylab("Het/Hom") + ggtitle("Het/Hom - INDELSs, WG") + geom_text(aes(label=IDhet_hom)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))
p4a <- ggplot(aes(y = nHets, x=1:nrow(newdata4o)), data=newdata4o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. heterozygous") + theme( plot.margin=unit(c(0.5,0.5,-1,0), "cm"))
p4b <- ggplot(aes(y = nHom, x=1:nrow(newdata4o)), data=newdata4o) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of variants") + ylab("N. homozygous") + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))



grid.arrange(p1, p2, p3,p4, p1a, p2a, p3a,p4a, p1b, p2b, p3b,p4b, heights=c(2/16, 1/16, 1/16, 2/16, 1/16, 1/16,2/16,1/16,1/16,2/16,1/16,1/16), ncol=4)

dev.off()


newdata1DP <- merge(newdata1,DPd,by.x="IDinsdel",by.y="ID")
newdata1DP <- newdata1DP[order(newdata1DP$mean), ]

newdata2DP <- merge(newdata2,DPd,by.x="IDinsdel",by.y="ID")
newdata2DP <- newdata2DP[order(newdata2DP$mean), ]

newdata3DP <- merge(newdata3,DPd,by.x="IDsing",by.y="ID")
newdata3DP <- newdata3DP[order(newdata3DP$mean), ]

newdata4DP <- merge(newdata4,DPd,by.x="IDsing",by.y="ID")
newdata4DP <- newdata4DP[order(newdata4DP$mean), ]



png(paste0(outD,"nVariantLoci_by_sample.png"), width=800, height=800, type="cairo")

p1 <- ggplot(aes(y = variant, x=1:nrow(newdata1DP)), data=newdata1DP) + geom_point(size=4,alpha=0.8,aes(color=mean)) +xlab("Samples ordered by DP") + ylab("N. Non-ref Variants") + ggtitle("N. Non-ref Variants - SNPs, EXOMES") + geom_text(aes(label=IDvariant)) + geom_hline(yintercept=mean(newdata1$variant),colour="yellow", linetype=1) + scale_colour_gradient(low = "blue", high="red")
p2 <- ggplot(aes(y = variant, x=1:nrow(newdata2DP)), data=newdata2DP) + geom_point(size=4,alpha=0.8,aes(color=mean)) +xlab("Samples ordered by DP") + ylab("N. Non-ref Variants") + ggtitle("N. Non-ref Variants - SNPs, WG") + geom_text(aes(label=IDvariant)) + geom_hline(yintercept=mean(newdata2$variant),colour="yellow", linetype=1) + scale_colour_gradient(low = "blue", high="red")
p3 <- ggplot(aes(y = variant, x=1:nrow(newdata3DP)), data=newdata3DP) + geom_point(size=4,alpha=0.8,aes(color=mean)) +xlab("Samples ordered by DP") + ylab("N.  Non-ref Variants") + ggtitle("N. Non-ref Variants - INDELSs, EXOMES") + geom_text(aes(label=IDvariant)) + geom_hline(yintercept=mean(newdata3$variant),colour="yellow", linetype=1) + scale_colour_gradient(low = "blue", high="red")
p4 <- ggplot(aes(y = variant, x=1:nrow(newdata4DP)), data=newdata4DP) + geom_point(size=4,alpha=0.8,aes(color=mean)) +xlab("Samples ordered by DP") + ylab("N. Non-ref Variants") + ggtitle("N. Non-ref Variants - INDELSs, WG") + geom_text(aes(label=IDvariant)) + geom_hline(yintercept=mean(newdata4$variant),colour="yellow", linetype=1) + scale_colour_gradient(low = "blue", high="red")

multiplot(p1, p2, p3, p4, cols = 1)

dev.off()


newdata3a <- newdata3[order(newdata3$nindel),]
newdata4a <- newdata4[order(newdata4$nindel),]

png(paste0(outD,"insertionDeletionRatio_by_sample.png"), width=1400, height=800, type="cairo")


p3 <- ggplot(aes(y = insdel, x=1:nrow(newdata3a)), data=newdata3a) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of indels") + ylab("Insertions/deletions") + ggtitle("Insertions/deletions - INDELSs, EXOMES") + geom_text(aes(label=IDinsdel)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))

p3a <- ggplot(aes(y = nDeletions, x=1:nrow(newdata3a)), data=newdata3a) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of indels") + ylab("N. deletion") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))

p3b <- ggplot(aes(y = nInsertions, x=1:nrow(newdata3a)), data=newdata3a) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of indels") + ylab("N. insertions")  + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))


p4 <- ggplot(aes(y = insdel, x=1:nrow(newdata4)), data=newdata4a) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of indels") + ylab("Insertions/deletions") + ggtitle("Insertions/deletions - INDELSs, WG") + geom_text(aes(label=IDinsdel)) + geom_text(aes(label=IDinsdel)) + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))

p4a <- ggplot(aes(y = nDeletions, x=1:nrow(newdata4a)), data=newdata4a) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of indels") + ylab("N. deletion") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.margin=unit(c(0.5,0.5,-1,0.5), "cm"))

p4b <- ggplot(aes(y = nInsertions, x=1:nrow(newdata4a)), data=newdata4a) + geom_point(size=4,alpha=0.8) +xlab("Samples ordered by N. of indels") + ylab("N. insertions")  + theme( plot.margin=unit(c(0.5,0.5,0.5,0), "cm"))


grid.arrange(p3, p4, p3a,p4a, p3b,p4b, heights=c(2/8, 1/8, 1/8, 2/8, 1/8, 1/8), ncol=2)

dev.off()


png(paste0(outD,"TITV_by_sample.png"), width=800, height=800, type="cairo")

p1 <- ggplot(aes(y = titv, x=1:nrow(newdata1)), data=newdata1) + geom_point(size=4,alpha=0.8) +xlab("Samples") + ylab("Ti/Tv") + ggtitle("Ti/Tv - SNPs, EXOMES") + geom_text(aes(label=IDtitv))
p2 <- ggplot(aes(y = titv, x=1:nrow(newdata2)), data=newdata2) + geom_point(size=4,alpha=0.8) +xlab("Samples") + ylab("Ti/Tv") + ggtitle("Ti/Tv - SNPs, WG") + geom_text(aes(label=IDtitv))

multiplot(p1, p2, cols = 1)

dev.off()

