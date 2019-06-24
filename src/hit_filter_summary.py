import argparse
import sys
import os
import re
from file_handling_lib import *


# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes mapping from filtering blast output")
parser.add_argument("blast_output", nargs="+", help="blast output for which to summarize mapping")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
parser.add_argument("--use_type", "-u", action="store_true", help="Should FASTQ type be associated with stats?")
parser.add_argument("--use_sample", "-s", action="store_true", help="Should summary stats be associated with the sample ID instead of the blast file name?")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

column_names = []
# If we are given a sample ID, we associate this summary with the sample ID
if args.use_sample:
    column_names.append("sample")

# Otherwise, we use the name of the FASTQ
else:
    column_names.append("blast_output")

# Add the column names for the summary stats to the list of column names
summary_stat_names = ["post_hit_filter_matched_reads", "post_hit_filter_matches", "post_hit_filter_average_e_value"]

# Create column names for the summary stats
for name in summary_stat_names:

    # If we are associating the summary stats with the FASTQ type, indicate that in each summary stat column header
    if args.use_type:
        fastq_type = re.search("^[^.]*\\.([^.]*)\\.", os.path.basename(args.fastq)).group(1)
        column_names.append("_".join([fastq_type, name]))

    # Otherwise, use the plain summary stat column header
    else:
        column_names.append(name)

output.write("\t".join(column_names) + os.linesep)

# Initialize trackers for parsing a file
matched_read_count = 0
matched_reads = set()
match_count = 0
total_e_value = 0

# Parse the blast output file
f = custom_read(args.blast_output)
for line in f:

    # Each line is a match
    match_count += 1

    # Now check if we've found a match for a read that didn't already have a match
    split_line = line.strip().split("\t")
    read_name = split_line[0]
    e_value = float(split_line[2])
    total_e_value += e_value

    # If we've found the first match for a read, count that as another read that had a match
    if read_name not in matched_reads:
        matched_read_count += 1
        matched_reads.add(read_name)

row_id = os.path.basename(args.blast_output)

# If we are using the sample ID, replace that as the row ID
if args.use_sample:
    row_id = re.search("^([^.])\\.", os.path.basename(args.blast_output)).group(1)

# Print the mapping summary for this blast output file
output.write("\t".join([row_id, str(matched_read_count), str(match_count), str(total_e_value/match_count)]) + os.linesep)
