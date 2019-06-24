import argparse
import sys
import os
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Using a UNIREF ID mapping file, creates a table mapping UNIREF gene IDs to ortholog groups.")
parser.add_argument("uniref_id_mapping", help="File from UNIREF mapping gene IDs to various annotations")
parser.add_argument("ortholog_type", help="The identifier for the type of ortholog to map genes to (e.g. KO, eggNOG)")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Open the output stream
# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Scan the UNIREF ID mapping file
with custom_read(args.uniref_id_mapping) as full_mapping:
    for line in full_mapping:

        # Split the line into its three components (UNIREF ID, mapping type, mapped value)
        mapping_info = line.strip().split()

        # If the mapping type matches the requested type, add this mapping to the output (ortholog first column, gene second column)
        if mapping_info[1] == args.ortholog_type:
            output.write("\t".join([mapping_info[2], mapping_info[0]]) + os.linesep)

output.close()
