#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 11/21/2017

import argparse,sys,os.path
from file_handling import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes mapping from blast output")
parser.add_argument("blast_output_files", nargs="+", help="blast output files for which to summarize mapping")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the mapping summary file
output.write("\t".join(["blast_output_file", "matched_reads", "matches", "average_e_value"]) + "\n")

# Parse each blast output file
for filename in args.blast_output_files:

    # Initialize trackers for parsing a file
    matched_read_count = 0
    matched_reads = set()
    match_count = 0
    total_e_value = 0

    # Parse the blast output file
    f = custom_read(filename)
    for line in f:

        # Each line is a match
        match_count += 1

        # Now check if we've found a match for a read that didn't already have a match
        split_line = line.strip().split("\t")
        read_name = split_line[0]
        e_value = float(split_line[10])
        total_e_value += e_value

        # If we've found the first match for a read, count that as another read that had a match
        if read_name not in matched_reads:
            matched_read_count += 1
            matched_reads.add(read_name)

    # Print the mapping summary for this blast output file
    output.write("\t".join([os.path.basename(filename), str(matched_read_count), str(match_count), str(total_e_value/match_count)]) + "\n")
