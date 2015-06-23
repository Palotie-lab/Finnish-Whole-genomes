#!/bin/bash
### PROGRAM TO CALCULATE ALLELE BALANCE AND % OF ALLELE DEVIATING FROM 30/70 ratio ###
### Argument 1: Input .txt file with AD and GT info
### Argument 2: Output folder where files should be saved

import numpy as np
import sys
from time import gmtime, strftime


def file_nsamp(fname):
    with open(fname) as f:
        for l in f:
            fields = l.split(" ")
            return int((len(fields)-1)/2)
        
def func(a, b):
  return not set(a).isdisjoint(b)


def file_len(fname):
    n_lines=0
    with open(fname) as f:
        for l in f:
            n_lines+=1
    return n_lines
    
    
het=["0/1","1/0"]


inD=sys.argv[1]
outD=sys.argv[2]

inlast = ("_").join(inD.split("_")[-2:])

#inD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/AD_G77318RH_table_SNP_EX.txt"
#outD="/humgen/atgu1/fs03/wip/aganna/fin_seq/processed/seq/temp/"

namesfile = inD

ABw=open(outD+'AB_'+ inlast, 'w') 
ABPw=open(outD+'ABP_'+ inlast, 'w')

f=open(namesfile)

n_samples=file_nsamp(namesfile)
n_lines=file_len(namesfile)

print '-' * 10
print 'The number of samples is %s' % n_samples
print '-' * 10

for q,line in enumerate(f):
    line = line.strip()
    fields = line.split(" ")
    assert(len(fields)==n_samples*2)
    if func(het,fields):
        AB=[]
        ABP=0
        for i in range(0, n_samples*2, 2):
            if fields[i] == "0/1" or fields[i] == "1/0":
                AD=fields[i+1].split(",")
                if (np.float64(AD[1])!=0):
                    ab=np.float64(AD[0])/(np.float64(AD[0])+np.float64(AD[1]))
                    ar=np.float64(AD[0])/np.float64(AD[1])
                    AB.append(ab)
                    if ar < (3/7) or ar > (7/3):
                        ABP+=1         
        if len(AB) != 0:
            towrite=np.float(ABP)/np.float(len(AB))
            ABPw.write(str(towrite))
            ABPw.write("\n")
            ABw.write(str(round(np.nanmean(np.array(AB)),2)))
            ABw.write("\n") 
        if (q % 1000000) == 0:
            print('Processing line',q,'out of',n_lines)
            print strftime("%Y-%m-%d %H:%M:%S", gmtime())
ABw.close()
ABPw.close()

