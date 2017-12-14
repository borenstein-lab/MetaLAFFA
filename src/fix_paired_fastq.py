#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: AlexE
# Date: 12/15/2014

NUM_LINES_PER_READ=4

import argparse
from file_handling import *
from future import *

parser = argparse.ArgumentParser()
parser.add_argument("fastq", help="One of the paired end fastq files to fix sequence names for.", metavar="Fastq")
parser.add_argument("tag", help="The tag to add to the end of sequence tags so reads are properly paired", metavar="Tag")
args = parser.parse_args()

f = custom_read(args.fastq)
line = f.readline()
count = 0
while line != "":
    if count == 0:
        split_line = line.strip().split()
        if re.search("/[0-9]+", split_line[0]):
            print(split_line[0][:-1] + args.tag)
        elif split_line[0][-1] == '/':
            print(split_line[0] + args.tag)
        else:
            print(split_line[0] + '/' + args.tag)
    else:
        print(line.strip())
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()
