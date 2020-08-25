#!/usr/bin/env python

import argparse
import sys
import os
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Using a UniRef ID mapping file, creates a table mapping UniRef gene IDs to ortholog groups.")
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

# Scan the UniRef ID mapping file
id_to_labels = {}
with custom_read(args.uniref_id_mapping) as full_mapping:
    for line in full_mapping:

        # Split the line into its three components (ID, mapping type, mapped value)
        mapping_info = line.strip().split()

        # If the mapping type matches the requested UniRef version, add this mapping to the id-to-uniref id mapping
        if mapping_info[1].lower() == args.uniref_version:
            if mapping_info[0] not in id_to_labels:
                id_to_labels[mapping_info[0]] = {"uniref": set(), "ortholog": set()}
            id_to_labels[mapping_info[0]]["uniref"].add(mapping_info[2])

        # If the mapping type matches the requested ortholog type, add this mapping to the id-to-ortholog mapping
        if mapping_info[1].lower() == args.ortholog_type:
            if mapping_info[0] not in id_to_labels:
                id_to_labels[mapping_info[0]] = {"uniref": set(), "ortholog": set()}
            id_to_labels[mapping_info[0]]["ortholog"].add(mapping_info[2])

# Convert mapping from IDs to UniRef ID and ortholog to an ortholog-to-UniRef ID mapping
ortholog_to_uniref = {}
ids = list(id_to_labels.keys())
for curr_id in ids:
    if len(id_to_labels[curr_id]["uniref"]) > 0 and len(id_to_labels[curr_id]["ortholog"]) > 0:
        for ortholog in id_to_labels[curr_id]["ortholog"]:
            if ortholog not in ortholog_to_uniref:
                ortholog_to_uniref[ortholog] = set()
            for uniref in id_to_labels[curr_id]["uniref"]:
                ortholog_to_uniref[ortholog].add(uniref)

    # Remove the entry from the first map to save memory
    del id_to_labels[curr_id]

# For each ortholog, write a line for each UniRef version ID that maps to it
for ortholog in ortholog_to_uniref.keys():
    for uniref in ortholog_to_uniref[ortholog]:
        output.write("\t".join([ortholog, uniref]) + os.linesep)

output.close()
