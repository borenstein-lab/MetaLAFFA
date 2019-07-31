import argparse
import sys
from file_handling_lib import *

NUM_LINES = 4

parser = argparse.ArgumentParser()
parser.add_argument("fastq1", help="The fastq file with the first half of the paired reads.", metavar="Fastq 1")
parser.add_argument("fastq2", help="The fastq file with read names that need to be altered.", metavar="Fastq 2")
args = parser.parse_args()

f1 = custom_read(args.fastq1)

f2 = custom_read(args.fastq2)

line1 = f1.readline()
line2 = f2.readline()
count = 0
while line1 != "" and line2 != "":
    if count == 0 or count == 2:
        print(line1.strip())
    else:
        print(line2.strip())
    count += 1
    count = count % 4
    line1 = f1.readline()
    line2 = f2.readline()
if line1 != "" or line2 != "":
    sys.stderr.write("WARNING: Numbers of reads did not match between files.  Reads may not be paired properly.\n")
f1.close()
f2.close()
