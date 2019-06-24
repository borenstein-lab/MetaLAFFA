import argparse
import os
import sys
import re
from fastq_summary_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Generates summary statistics about a quality-filtered FASTQ file")

# Add positional arguments
parser.add_argument("fastq", help="Fastq file for which to generate summary statistics")

# Add optional arguments
parser.add_argument("--output", "-o", help="Where output should be written (default: print to standard output)")
parser.add_argument("--use_type", "-u", action="store_true", help="Should FASTQ type be associated with stats?")
parser.add_argument("--use_sample", "-s", action="store_true", help="Should summary stats be associated with the sample ID instead of the FASTQ file name?")

args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Generate and the column names for the host filtering summary file
column_names = []

# If we are given a sample ID, we associate this summary with the sample ID
if args.use_sample:
    column_names.append("sample")

# Otherwise, we use the name of the FASTQ
else:
    column_names.append("fastq")

# Add the column names for the summary stats to the list of column names
summary_stat_names = ["post_quality_filter_reads", "post_quality_filter_average_read_length", "post_quality_filter_average_base_quality"]

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

# Get the stats for the FASTQ file
read_count, average_read_length, average_base_quality = get_fastq_stats(args.fastq)

row_id = os.path.basename(args.fastq)

# If we are using the sample ID, replace that as the row ID
if args.use_sample:
    row_id = re.search("^([^.])\\.", os.path.basename(args.fastq)).group(1)

# Print the host filtering summary for this file
output.write("\t".join([row_id, str(read_count), str(average_read_length), str(average_base_quality)]) + os.linesep)

output.close()
