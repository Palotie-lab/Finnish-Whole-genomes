### PLOT STAT CONCORDANCE CHIP-WGS ###
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

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/plots/"
# Read data


d <- read.table("/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/psych/concordance_SNP_EX", blank.lines.skip = F, skip=3, header=T, stringsAsFactor=F)

start <- which(d$Sample=="Sample")[1]
end <- which(d$Sample=="Sample")[2]

dsel <- d[(start+1):(end-4),]
dsel2 <- dsel[dsel$Sample!="ALL" & !dsel$Eval_Genotype%in%c("UNAVAILABLE","MIXED","Mismatching_Alleles") & !dsel$Comp_Genotype%in%c("UNAVAILABLE","MIXED","Mismatching_Alleles"),]

ConcHET <- NULL
ConcHOMVAR <- NULL
Non_ref_sens <- NULL
Non_ref_disc <- NULL
for (i in unique(dsel2$Sample))
{
  dseT <- dsel2[dsel2$Sample==i,]
  den <- sum(as.numeric(dseT$Proportion))
  concHET <- as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HET" & dseT$Eval_Genotype=="HET"]) / sum(as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HET"]))

  concHOMVAR <- as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HOM_VAR" & dseT$Eval_Genotype=="HOM_VAR"]) / sum(as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HOM_VAR"]))

  non_ref_sens <- (as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HOM_VAR" & dseT$Eval_Genotype=="HOM_VAR"]) +  as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HET" & dseT$Eval_Genotype=="HET"])) / sum(as.numeric(dseT$Proportion[dseT$Comp_Genotype=="HOM_VAR" | dseT$Comp_Genotype=="HET"]))

  non_ref_disc <-  sum(as.numeric(dseT$Proportion[!(dseT$Comp_Genotype=="HOM_REF" & dseT$Eval_Genotype=="HOM_REF") & !(dseT$Comp_Genotype=="NO_CALL" | dseT$Eval_Genotype=="NO_CALL") & !(dseT$Comp_Genotype=="HOM_VAR" & dseT$Eval_Genotype=="HOM_VAR") & !(dseT$Comp_Genotype=="HET" & dseT$Eval_Genotype=="HET")])) /
  sum(as.numeric(dseT$Proportion[!(dseT$Comp_Genotype=="HOM_REF" & dseT$Eval_Genotype=="HOM_REF") & !(dseT$Comp_Genotype=="NO_CALL" | dseT$Eval_Genotype=="NO_CALL")]))

  ConcHET <- c(ConcHET,concHET)
  ConcHOMVAR <- c(ConcHOMVAR,concHOMVAR)
  Non_ref_sens <- c(Non_ref_sens,non_ref_sens)
  Non_ref_disc <- c(Non_ref_disc,non_ref_disc)
}




newdata <- data.frame(ConcHET=ConcHET,ConcHOMVAR=ConcHOMVAR,Non_ref_sens=Non_ref_sens,Non_ref_disc=Non_ref_disc)

pdf("")


