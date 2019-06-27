import argparse
import sys
import os
import re
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes ortholog counting")
parser.add_argument("ortholog_profile", help="Ortholog profile for which to summarize ortholog counting")
parser.add_argument("--ortholog_name", "-n", default="ortholog", help="The name of the ortholog system used (e.g. KO, COG, etc.), will be used as the column header for ortholog IDs (default: %(default)s)")
parser.add_argument("--use_sample", "-s", action="store_true", help="Should summary stats be associated with the sample ID instead of the blast file name?")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
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
    column_names.append("ortholog_profile_file")

column_names.append(args.ortholog_name + "s")

# Print the column names for the gene counting summary file
output.write("\t".join(column_names) + os.linesep)

# Initialize trackers for parsing a file
ortholog_count = 0

# Parse the ortholog profile file
f = custom_read(args.ortholog_profile)

# Skip the column headers
line = f.readline()

# Just count all the lines of the table (different orthologs)
for line in f:
    ortholog_count += 1

row_id = os.path.basename(args.ortholog_profile)

# If we are using the sample ID, replace that as the row ID
if args.use_sample:
    row_id = re.search("^([^.]*)\\.", os.path.basename(args.ortholog_profile)).group(1)

# Print the KO counting summary for this file
output.write("\t".join([row_id, str(ortholog_count)]) + os.linesep)

output.close()
