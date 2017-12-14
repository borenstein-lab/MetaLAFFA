#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 12/14/2017

NUM_LINES_PER_READ = 4

import argparse,sys,os.path
from fastq_summarizing import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes the input read file")
parser.add_argument("fastqs", nargs="+", help="Fastq files for which to summarize read information")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the host filtering summary file
output.write("\t".join(["fastq_file", "reads", "average_read_length", "average_base_quality"]) + "\n")

# Parse each filtered fastq
for filename in args.fastqs:

    # Get the stats for the FASTQ file
    read_count, average_read_length, average_base_quality = get_fastq_stats(filename)

    # Print the host filtering summary for this file
    output.write("\t".join([os.path.basename(filename), str(read_count), str(average_read_length), str(average_base_quality)]) + "\n")
output.close()
