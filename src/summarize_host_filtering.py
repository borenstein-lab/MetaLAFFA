#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 11/30/2017

NUM_LINES_PER_READ = 4

import argparse,sys
from file_handling import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes host read filtering")
parser.add_argument("filtered_fastqs", nargs="+", help="Filtered fastq files for which to summarize host filtering")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the host filtering summary file
output.write("\t".join(["filtered_fastq_file", "post_host_filtering_reads"]) + "\n")

# Parse each filtered fastq
for filename in args.filtered_fastqs:

    # Initialize trackers for parsing a file
    read_count = 0

    # Parse the filtered fastq file
    f = custom_read(filename)
    count = 0
    line = f.readline()
    while line != "":

        # If we're on the first line of a read, then we found another read
        if count == 0:
            read_count += 1

        # Now update our counter for current line in a read's info
        count += 1
        count = count % NUM_LINES_PER_READ
        line = f.readline()

    # Print the host filtering summary for this file
    output.write("\t".join([filename, str(read_count)]) + "\n")
output.close()