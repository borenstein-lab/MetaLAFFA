#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: AlexE
# Date: 12/30/2014

import argparse, re, gzip

parser = argparse.ArgumentParser()
parser.add_argument("marked", help="The file of marked reads to remove.", metavar="Marked")
parser.add_argument("fastq", help="The fastq file to remove marked reads from.", metavar="Fastq")
args = parser.parse_args()

# Read in the marked reads
marked_file = None
if re.search("\.gz$", args.marked):
    marked_file = gzip.open(args.marked, 'r')
else:
    marked_file = open(args.marked, 'r')
marked_reads = set()
for line in marked_file:
    read_name = line.strip()
    if read_name[-2] == '/':
        marked_reads.add(read_name[:-2])
    else:
        marked_reads.add(read_name)
marked_file.close()

fastq = None
if re.search("\.gz$", args.fastq):
    fastq = gzip.open(args.fastq, 'r')
else:
    fastq = open(args.fastq, 'r')
line = fastq.readline()
while line != "":
    read_name = line.strip().split()[0][1:]
    if read_name[-2] == '/':
        read_name = read_name[:-2]
    if read_name in marked_reads:
        line = fastq.readline()
        line = fastq.readline()
        line = fastq.readline()
    else:
        print line.strip()
        line = fastq.readline()
        print line.strip()
        line = fastq.readline()
        print line.strip()
        line = fastq.readline()
        print line.strip()
    line = fastq.readline()
fastq.close()
