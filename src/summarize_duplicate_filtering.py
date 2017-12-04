#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 11/30/2017

NUM_LINES_PER_READ = 4

import argparse,sys
from fastq_summarizing import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes duplicate read filtering")
parser.add_argument("filtered_fastqs", nargs="+", help="Filtered fastq files for which to summarize duplicate filtering")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the duplicate filtering summary file
output.write("\t".join(["filtered_fastq_file", "post_duplicate_filtering_reads", "post_duplicate_filtering_average_read_length", "post_duplicate_filtering_average_base_quality"]) + "\n")

# Parse each filtered fastq
for filename in args.filtered_fastqs:

    # Get the stats for the FASTQ file
    read_count, average_read_length, average_base_quality = get_fastq_stats(filename)

    # Print the duplicate filtering summary for this file
    output.write("\t".join([filename, str(read_count), str(average_read_length), str(average_base_quality)]) + "\n")
output.close()