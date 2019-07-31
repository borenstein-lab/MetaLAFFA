import argparse
import re
import gzip
from file_handling_lib import *

parser = argparse.ArgumentParser()
parser.add_argument("marked", help="The file of marked reads to remove.", metavar="Marked")
parser.add_argument("fastq", help="The fastq file to remove marked reads from.", metavar="Fastq")
args = parser.parse_args()

# Read in the marked reads
marked_file = custom_read(args.marked)
marked_reads = set()
for line in marked_file:
    read_name = line.strip()
    # Remove the part of the read name identifying which end of a read pair it is from if present
    end_identifier_match = re.search("/[^/]*$", read_name)
    if end_identifier_match:
        read_name = read_name[:end_identifier_match.start()]
    marked_reads.add(read_name)
marked_file.close()

fastq = custom_read(args.fastq)
line = fastq.readline()
while line != "":
    read_name = line.strip().split()[0][1:]
    # Remove the part of the read name identifying which end of a read pair it is from if present
    end_identifier_match = re.search("/[^/]*$", read_name)
    if end_identifier_match:
        read_name = read_name[:end_identifier_match.start()]
    if read_name in marked_reads:
        line = fastq.readline()
        line = fastq.readline()
        line = fastq.readline()
    else:
        print(line.strip())
        line = fastq.readline()
        print(line.strip())
        line = fastq.readline()
        print(line.strip())
        line = fastq.readline()
        print(line.strip())
    line = fastq.readline()
fastq.close()
