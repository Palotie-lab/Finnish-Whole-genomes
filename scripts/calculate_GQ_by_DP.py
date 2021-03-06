#!/bin/bash
### PLOT DP BY GQ GENOTYPE-SPECIFIC METRICS ###
### Argument 1: a file generated from bcftool including DP genotype-specific values 
###	e.g. /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/DP_G77318RH_table_SNP_EX.txt
### Argument 2: a file generated from bcftool including GQ genotype-specific values 
###	e.g. /humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/GQ_G77318RH_table_SNP_EX.txt
### Argument 3: a path were to ouput the files


import sys
import pandas as pd
import numpy as np
from time import gmtime, strftime


in1=sys.argv[1]
in2=sys.argv[2]
in3=sys.argv[3]
out=sys.argv[4]

#in1='/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/DP_G89387_INDEL_EX.txt'
#in2='/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/GQ_G77318RH_INDEL_EX.txt'
#in3='/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/samplenames_G77318RH_INDEL_EX.txt'
#out='/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/'

inlast = ("_").join(in1.split("_")[-3:])


print "Loading DP values"
dp1 = pd.read_csv(in1 ,delim_whitespace=True,engine='c', header=None, na_values='.')

print "Loading GQ values"
gq1 = pd.read_csv(in2 ,delim_whitespace=True,engine='c', header=None, na_values='.')

names1= pd.read_csv(in3 ,delim_whitespace=True,engine='c', header=None, na_values='.')


def f(x):
	return float((x[np.isfinite(x)]>30).sum())/float(len(x[np.isfinite(x)]))


def GQGT20byDP(gq,dp):
    Gqgt20dplt10=[]
    Dagggt=[]
    Dagglt=[]
    for i in xrange(len(gq.columns)):
        gqt=gq.ix[:,i]
        dpt=dp.ix[:,i]
        gqtnm = gqt[np.isfinite(gqt) & np.isfinite(dpt)]
        dptnm = dpt[np.isfinite(gqt) & np.isfinite(dpt)]
        assert len(gqtnm)==len(dptnm),"ERROR: different length"
        tf = gqtnm > 20  
        #gqtnm20 = gqtnm < 20
        #dptnm30 = dptnm > 30
        #gqgt20dplt10 = (gqtnm20 & dptnm30).sum() 
        agggt=[]
        agglt=[]
        if (i % 100) == 0:
    	    print('Processing sample',i,'out of',len(gq1.loc[1,]))
    	    print strftime("%Y-%m-%d %H:%M:%S", gmtime())
        for j in xrange(5,80,5):
    	    dflt = dptnm < j
            dfgt = (dptnm < j) & (dptnm > j-5)
    	    newnlt = (tf & dflt).sum()
            newngt = (tf & dfgt).sum()
            if dfgt.sum() == 0:
                agglt.append('NA')
                agggt.append('NA')
            else:
                agglt.append(float(newnlt)/float(len(tf)))
                agggt.append(float(newngt)/float(dfgt.sum()))
                          
        Dagggt.append(agggt)
        Dagglt.append(agglt)
        #Gqgt20dplt10.append(gqgt20dplt10)
    return(Dagggt,Dagglt)

print "Processing DGQ > 20 by DP combinations"
Dagggt, Dagglt = GQGT20byDP(gq1,dp1)

GQgt20byDPlt = pd.DataFrame(Dagglt).transpose()
GQgt20byDPgt = pd.DataFrame(Dagggt).transpose()

#Gqgt20dplt10F = pd.concat([names1[[0]],pd.DataFrame(Gqgt20dplt10)],axis=1)

dpgt30=dp1.apply(f,axis=0)
GQmean=pd.concat([names1[[0]],gq1.mean(axis=0)],axis=1)

dpgt30.to_csv(out + 'DPgt30_' + inlast,sep=' ')

GQgt20byDPlt.to_csv(out + 'GQgt20byDPlt_' + inlast,sep=' ', header=[l[0] for l in names1[[0]].values.tolist()])
GQgt20byDPgt.to_csv(out + 'GQgt20byDPgt_' + inlast,sep=' ', header=[l[0] for l in names1[[0]].values.tolist()])

#Gqgt20dplt10F.to_csv(out + 'Gqgt20dplt10_' + inlast,sep=' ',header=None)
GQmean.to_csv(out + 'GQmean_' + inlast,sep=' ',header=None)

