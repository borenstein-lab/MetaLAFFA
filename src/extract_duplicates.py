#!/usr/bin/env python

import argparse

NAME_FIELD = 0
FLAG_FIELD = 1
DUP_FLAG = -11

parser = argparse.ArgumentParser()
parser.add_argument("samview", help="The samview of a SAM or BAM file with marked duplicate reads.", metavar="Samview")
parser.add_argument("-s", "--suffix", help="The suffix to add to read names if necessary (for paired reads, for example).", metavar="Suffix")
args = parser.parse_args()

f = open(args.samview, 'r')
line = f.readline()
while line != "":
    fields = line.strip().split('\t')

    # Grab the duplicate flag
    flag = bin(int(fields[FLAG_FIELD]))
    if len(flag) >= 13:
        if flag[DUP_FLAG] == "1":
            if args.suffix:
                print(fields[NAME_FIELD] + args.suffix)
            else:
                print(fields[NAME_FIELD])

    # Skip the next line, which is the paired read, and move to the next read
    line = f.readline()
    line = f.readline()
f.close()
