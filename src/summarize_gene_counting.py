#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 11/30/2017

import argparse,sys,os.path
from file_handling import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes gene counting")
parser.add_argument("gene_profile_files", nargs="+", help="Gene profile files for which to summarize gene counting")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the gene counting summary file
output.write("\t".join(["gene_profile_file", "genes"]) + "\n")

# Parse each gene profile
for filename in args.gene_profile_files:

    # Initialize trackers for parsing a file
    gene_count = 0

    # Parse the gene profile file
    f = custom_read(filename)

    # Skip the column headers
    line = f.readline()

    # Just count all the lines of the table (different genes)
    for line in f:
        gene_count += 1

    # Print the gene counting summary for this file
    output.write("\t".join([os.path.basename(filename), str(gene_count)]) + "\n")
output.close()