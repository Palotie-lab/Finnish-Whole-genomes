__author__ = 'aganna'

import gzip
import argparse
import numpy as np
import os
import logging
import time
import sys
#import line_profiler

def file_len(file):
    with gzip.open(file) if args.vcf.endswith('.gz') else open(file) as f:
        data_lines = 0
        for line in f:
            if line[0] != "#":
                data_lines += 1
    return data_lines

class QuickTimer():
    def __init__(self):
        self.timer = time.time()

    def print_time(self):
        return str(round(time.time() - self.timer, ndigits=2))

    def get_time(self):
        return time.time() - self.timer

    def reset(self):
        self.timer = time.time()
  


# def get_ID(file,LEN_SITES):
#     IDF = np.empty(LEN_SITES, dtype=np.dtype('a100'))
#     with gzip.open(file, 'r') as f:
#         for num,line in enumerate(f):
#             if line[0] != "#":
#                 fields = line.strip().split('\t')
#                 ID = fields[0] + ":" + fields[1] + ":" + fields[3] + ":" + fields[4]
#                 IDF[num]=ID
#     return IDF


#@profile
def main(args):
    logging_format = "%(levelname)s: %(message)s"
    logging.basicConfig(level=logging.DEBUG, format=logging_format)

    f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)
    g = open(args.output, 'w')
    e = gzip.open(os.path.splitext(args.output)[0] + "_exceptions.txt.gz", 'w')

    print "getting length"
    LEN_SITES = file_len(args.vcf)

    # print "get IDF"
    # IDF = get_ID(args.vcf,LEN_SITES)
    # print IDF[0:10]

    logging.info("Starting site processing")
    timer = QuickTimer()
    prog = -1

    for num, line in enumerate(f):

        if ((num+1) * 100) / LEN_SITES > prog:
            prog += 1
            if prog > 1:
                sys.stdout.write("\r\033[KPercent complete: %d%%  Elapsed time: %s seconds / %d minutes remaining" %
                                 (prog, timer.print_time(), int(float(timer.get_time()) * (100-prog)/prog / 60)))
                sys.stdout.flush()
            else:
                sys.stdout.write("\rPercent complete: %d%%" % prog)
                sys.stdout.flush()


        if line.startswith('##'):
            g.write(line)
            #continue
        elif line.startswith('#'):
            lineT = line.split()
            fmt_index = lineT.index("FORMAT") 
            nsamp = len(lineT[fmt_index+1:])
            maxsamp = 0.20*nsamp
            print "Working with " + str(nsamp) + " samples"
            g.write(line)
        else:
            fields = line.strip().split('\t')
            ID = fields[0] + ":" + fields[1] + ":" + fields[3] + ":" + fields[4]

            fmt = fields[fmt_index]
            gt_index = fmt.split(":").index("GT")
            ad_index = fmt.split(":").index("AD")
            dp_index = fmt.split(":").index("DP")
            gq_index = fmt.split(":").index("GQ")

            alt = fields[4]
            n_alleles = len(alt.split(",")) + 1
            gt_countM = 0
            gt_countGQ = 0
            gt_countAD = 0
            gt_countDP = 0

            newfield = [None] * nsamp

            for ind, sample_field in enumerate(fields[fmt_index+1:]):
                sample_split = sample_field.split(":")

                gt = sample_split[gt_index]

                if gt == "./.":
                    gt_countM +=1
                    newfield[ind] = sample_field
                    continue

                ad = sample_split[ad_index]
                dp = sample_split[dp_index]
                if dp == ".":
                    dp = 700
                gq = sample_split[gq_index]

                all_ba = 0.5
                if  gt.split("/")[0] != gt.split("/")[1]:

                    ads = map(int, sample_split[ad_index].split(","))
                    ad_ref = ads[0]

                    if n_alleles == 2:
                        ad_alt = ads[-1]
                    else:
                        ad_alt = np.sum(ads[1:-1])

                    if (ad_alt+ad_ref) > 10:
                        all_ba = float(ad_alt) / float(ad_ref + ad_alt)

                if int(gq) < 20:
                    gt_countGQ +=1
                    gt = "./."

                if int(dp) > 400:
                    gt_countDP +=1
                    gt = "./."

                if all_ba < 0.20 or all_ba > 0.80:
                    gt_countAD +=1
                    gt = "./."

                sample_split[gt_index] = gt
                #newfield[ind] = str(gt) + ":" + str(ad) + ":" + str(dp) + ":" + str(gq) + ":" + str(pl)
                newfield[ind] = ':'.join(sample_split)
                
            gt_count = [gt_countM, gt_countGQ, gt_countDP, gt_countAD]

            if (sum(gt_count)) > maxsamp:
                e.write(str(ID) + " " + ' '.join(str(x) for x in gt_count) + '\n')
                continue

            assert len(newfield) == nsamp

            g.write('\t'.join(fields[0:fmt_index+1]) + '\t' + '\t'.join(x for x in newfield) + '\n')

    f.close()
    g.close()
    e.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('vcf',  type=str, help='Input VCF file; may be gzipped')
    #parser.add_argument('hwe',  type=str, help='File from plink containing the HWE results')
    parser.add_argument('output', type=str, help='Output file, may be gzipped if ends in .gz')
    args = parser.parse_args()
    main(args)
