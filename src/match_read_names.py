#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: AlexE
# Date: 12/5/2014

import argparse, re, sys, gzip

NUM_LINES = 4

parser = argparse.ArgumentParser()
parser.add_argument("fastq1", help="The fastq file with the first half of the paired reads.", metavar="Fastq 1")
parser.add_argument("fastq2", help="The fastq file with read names that need to be altered.", metavar="Fastq 2")
args = parser.parse_args()

# Adrian: replace with your file_handling functions
f1 = None
if re.search("\.gz$", args.fastq1):
    f1 = gzip.open(args.fastq1, 'r')
else:
    f1 = open(args.fastq1, 'r')

# Adrian: replace with your file_handling functions
f2 = None
if re.search("\.gz$", args.fastq2):
    f2 = gzip.open(args.fastq2, 'r')
else:
    f2 = open(args.fastq2, 'r')

line1 = f1.readline()
line2 = f2.readline()
count = 0
while line1 != "" and line2 != "":
    if count == 0 or count == 2:
        print line1.strip()
    else:
        print line2.strip()
    count += 1
    count = count % 4
    line1 = f1.readline()
    line2 = f2.readline()
if line1 != "" or line2 != "":
    sys.stderr.write("WARNING: Numbers of reads did not match between files.  Reads may not be paired properly.\n")
f1.close()
f2.close()
