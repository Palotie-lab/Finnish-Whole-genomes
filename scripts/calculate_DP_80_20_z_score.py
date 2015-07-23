__author__ = 'aganna'
import logging
import argparse
import pandas
import gzip
import numpy as np
import math
import time
import sys
import random

class QuickTimer():
    def __init__(self):
        self.timer = time.time()

    def print_time(self):
        return str(round(time.time() - self.timer, ndigits=2))

    def get_time(self):
        return time.time() - self.timer

    def reset(self):
        self.timer = time.time()
  
def file_len(file):
    with gzip.open(file, 'r') as f:
        data_lines = 0
        for line in f:
            if line[0] != "#":
                data_lines += 1
    return data_lines
    
            
def main(args):
    logging_format = "%(levelname)s: %(message)s"
    logging.basicConfig(level=logging.DEBUG, format=logging_format)
    timer = QuickTimer()

    # Number of sites
    N_SITES = file_len(args.file)
    
    if  args.samplines > N_SITES:
        print 'Number of lines to sample > N. of lines: using maximum number of lines'
        samplines = N_SITES
    else:
        samplines = args.samplines

    
    with gzip.open(args.file, 'r') as f:
        

        header_line = "##"
        while header_line.startswith("##"):
            header_line = f.readline()
        
        
        # Select random lines #
        random_linenos = sorted(random.sample(xrange(N_SITES), samplines), reverse = True)
        
        lineno = random_linenos.pop()
        
         
        fmt_index = header_line.split().index("FORMAT")    
        sampnames = header_line.split("\t")[fmt_index+1:]
        num_samples = len(sampnames)
        all_bal = np.zeros(num_samples, dtype=float)
        abcollect = np.zeros(num_samples + 2, dtype=float)
        all_80_20 = np.zeros(num_samples)
        all_196 = np.zeros(num_samples)
        denominator= np.zeros(num_samples)
        zcollect=np.zeros(num_samples + 2, dtype=float)
        ref_reads = np.zeros(num_samples, dtype=np.int16)
        alt_reads = np.zeros(num_samples, dtype=np.int16)
         
        
        logging.info("Starting site processing")
        timer = QuickTimer()
        prog = -1
        
        
        
        for num, line in enumerate(f):
            
            
            # If is one of the randomly selected line
            if num == lineno:
            
                # Set vector for z scores to missing
                zcollect[:] = np.nan
                abcollect[:] = np.nan
                
                # Set up the timer
                if ((num+1) * 100) / N_SITES > prog:
                    prog += 1
                    if prog > 1:
                        sys.stdout.write("\r\033[KPercent complete: %d%%  Elapsed time: %s seconds / %d minutes remaining" %
                                         (prog, timer.print_time(), int(float(timer.get_time()) * (100-prog)/prog / 60)))
                        sys.stdout.flush()
                    else:
                        sys.stdout.write("\rPercent complete: %d%%" % prog)
                        sys.stdout.flush()
                    
                
                # Split line
                line = line.strip().split("\t")
                
                # When exausting number of lines than stop the pgm
                if len(random_linenos) > 0:
                    lineno = random_linenos.pop()
                else:
                    break
            
                # X chromosome == 23
                if line[0] == "X":
                    zcollect[0] = 23
                    abcollect[0] = 23
                elif line[0] == "Y":
                    zcollect[0] = 24
                    abcollect[0] = 24
                elif line[0] == "MT":
                    zcollect[0] = 25
                    abcollect[0] = 25
                else:
                    zcollect[0] = line[0]
                    abcollect[0] = line[0]
            
                zcollect[1] = line[1]
                abcollect[1] = line[1]
                
                        
                ref, alt = line[3], line[4]
                n_alleles_line = len(alt.split(",")) + 1

            
                fmt = line[fmt_index]
                gt_index = fmt.split(":").index("GT")
                ad_index = fmt.split(":").index("AD")
            
                # Per-column calculations #
                for ind, sample_field in enumerate(line[fmt_index+1:]):
                  
                    sample_split = sample_field.split(":")
                    gt = sample_split[gt_index]

                    if gt == "./.":
                        continue
                        
                    ad = map(int, sample_split[ad_index].split(","))
                    ref_reads[ind] = ad[0]

                    if n_alleles_line == 2:
                        alt_reads[ind] = ad[-1]
                    else:
                        alt_reads[ind] = np.sum(ad[1:-1])

                    # Only heterozygous
                    if  gt.split("/")[0] != gt.split("/")[1] and gt.split("/")[0]=="0":
                        try:
                        
                            altad=alt_reads[ind]
                            refad=ref_reads[ind]
                            all_ba = float(altad) / float(refad + altad)
                                                    
                        except ZeroDivisionError:
                            continue
                        
                        # Z-scores for normal aprximation to bionomial dist
                        m=(refad+altad)*0.47
                        s=math.sqrt((refad+altad)*0.47*(1-0.47))
                        z=(altad-m)/s
                        
                        # Only if DP > 10 #
                        if m < 5:
                            continue
                        zcollect[ind+2]=z
                        abcollect[ind+2]=all_ba
                        all_bal[ind] += all_ba
                        denominator[ind] += 1  
                        
                        if z < (float(-1.96)) or z > (float(1.96)):
                            all_196[ind] += 1
                    
                toout=' '.join(map(str, np.round(zcollect,2)))
                toout2=' '.join(map(str, np.round(abcollect,2)))
                                 
                zwrite.write(toout)
                zwrite.write('\n')
            
                albalwrite.write(toout2)
                albalwrite.write('\n')
    
    print
    
    # Summary stats per sample #

    ab = pandas.DataFrame(sampnames)
    ab.columns = ['ID']
    ab['AB'] = all_bal/denominator
    ab['N_HET_DPgt10'] = denominator
    ab['AB_Z_196'] = all_196/denominator

    ab.to_csv(args.out + "_AB_stats_per_sample.txt", sep=",",index_label=False,index=False)
  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='calculate_DP_80_20_z_score.py',
                                     formatter_class=lambda prog:
                                     argparse.HelpFormatter(prog, max_help_position=30))

    parser.add_argument('file',
                        type=str,
                        help='- Path to VCF subset file')
    parser.add_argument('out',
                        type=str,
                        help='- Path + file suffix')
                        
    parser.add_argument('samplines',
                            type=int,
                            help='- Number of lines to sample')

    
    args, pass_through_args = parser.parse_known_args()
    zwrite = open(args.out + "_Zscores_per_sample.txt", 'w')
    albalwrite = open(args.out + "_Albal_per_sample.txt", 'w')
    main(args)
    zwrite.close()
    albalwrite.close()
    

