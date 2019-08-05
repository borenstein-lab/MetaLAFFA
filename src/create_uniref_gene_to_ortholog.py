import argparse
import sys
import os
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Using a UniRef ID mapping file, creates a table mapping UniRef gene IDs to ortholog groups.")
parser.add_argument("uniref_id_mapping", help="File from UniRef mapping gene IDs to various annotations")
parser.add_argument("uniref_version", help="The version of the UniRef database to match orthologs to (e.g. uniref90)")
parser.add_argument("ortholog_type", help="The identifier for the type of ortholog to map genes to (e.g. ko, eggnog)")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Standardize strings to be lower-case
args.uniref_version = args.uniref_version.lower()
args.ortholog_type = args.ortholog_type.lower()

# Open the output stream
# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Scan the UNIREF ID mapping file
id_to_uniref_version = {}
ortholog_to_id = {}
with custom_read(args.uniref_id_mapping) as full_mapping:
    for line in full_mapping:

        # Split the line into its three components (ID, mapping type, mapped value)
        mapping_info = line.strip().split()

        # If the mapping type matches the requested UniRef version, add this mapping to the id-to-uniref version mapping
        if mapping_info[1].lower() == args.uniref_version:
            if mapping_info[0] not in id_to_uniref_version:
                id_to_uniref_version[mapping_info[0]] = []
            id_to_uniref_version[mapping_info[0]].append(mapping_info[2])

        # If the mapping type matches the requested ortholog type, add this mapping to the ortholog-to-id mapping
        if mapping_info[1].lower() == args.ortholog_type:
            if mapping_info[2] not in ortholog_to_id:
                ortholog_to_id[mapping_info[2]] = []
            ortholog_to_id[mapping_info[2]].append(mapping_info[0])

# For each ortholog, write a line for each UniRef version ID that maps to it
for ortholog in ortholog_to_id.keys():
    for uniref_id in ortholog_to_id[ortholog]:
        if uniref_id in id_to_uniref_version:
            for uniref_version_id in id_to_uniref_version[uniref_id]:
                output.write("\t".join([ortholog, uniref_version_id]) + os.linesep)

output.close()
