#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: AlexE
# Date: 12/15/2014

NAME_FIELD = 0
FLAG_FIELD = 1
DUP_FLAG = -11

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("samview", help="The samview of a SAM or BAM file with marked duplicate reads.", metavar="Samview")
parser.add_argument("-s", "--suffix", help="The suffix to add to read names if necessary (for paired reads, for example).", metavar="Suffix")
args = parser.parse_args()

f = open(args.samview, 'r')
line = f.readline()
while line != "":
    fields = line.strip().split('\t')
    flag = bin(int(fields[FLAG_FIELD]))
    if len(flag) >= 13:
        if flag[DUP_FLAG] == "1":
            if args.suffix:
                print fields[NAME_FIELD] + args.suffix
            else:
                print fields[NAME_FIELD]
    # Adrian: Why are two lines being read in here?
    line = f.readline()
    line = f.readline()
f.close()
