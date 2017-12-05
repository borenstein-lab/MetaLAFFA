#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 11/30/2017

import argparse,sys
from file_handling import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes KO counting")
parser.add_argument("KO_profile_files", nargs="+", help="Gene profile files for which to summarize KO counting")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the KO counting summary file
output.write("\t".join(["KO_profile_file", "KOs"]) + "\n")

# Parse each KO profile
for filename in args.KO_profile_files:

    # Initialize trackers for parsing a file
    KO_count = 0

    # Parse the KO profile file
    f = custom_read(filename)

    # Skip the column headers
    line = f.readline()

    # Just count all the lines of the table (different KOs)
    for line in f:
        KO_count += 1

    # Print the KO counting summary for this file
    output.write("\t".join([filename, str(KO_count)]) + "\n")
output.close()
