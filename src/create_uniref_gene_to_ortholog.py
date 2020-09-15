#!/usr/bin/env python

import argparse
import sys
import os
import re
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

# Read in the mapping file, identifying IDs that have definitions for both the specified UniRef version and the
# specified ortholog type. We only grab these IDs first so that we know which IDs we care about and can save memory
# down the line when we read in the data we use to match up the mapping types. This costs us some time because we
# will be scanning the file twice.
uniref_ids = set()
ortholog_ids = set()
paired_ids = set()
with custom_read(args.uniref_id_mapping) as full_mapping:
    for line in full_mapping:

        # Split the line into its three components (ID, mapping type, mapped value)
        mapping_info = line.strip().split()

        # Parse out the base ID (ignore isoform specification) and grab the mapping type
        mapping_id = re.match("^([^-]*)", mapping_info[0]).groups()[0]
        mapping_type = mapping_info[1].lower()

        # If this ID has not been identified as matching both desired mapping types, check for matches
        if mapping_id not in paired_ids:
            if mapping_type == args.uniref_version:
                if mapping_id in ortholog_ids:
                    paired_ids.add(mapping_id)
                    ortholog_ids.remove(mapping_id)
                else:
                    uniref_ids.add(mapping_id)
            elif mapping_type == args.ortholog_type:
                if mapping_id in uniref_ids:
                    paired_ids.add(mapping_id)
                    uniref_ids.remove(mapping_id)
                else:
                    ortholog_ids.add(mapping_id)

# Remove sets of IDs that don't match both mapping types to save memory
del uniref_ids
del ortholog_ids

# Now rescan the mapping file, extracting matching mapping type data for each ID
id_to_labels = {}
with custom_read(args.uniref_id_mapping) as full_mapping:
    for line in full_mapping:

        # Split the line into its three components (ID, mapping type, mapped value)
        mapping_info = line.strip().split()

        # Parse out the base ID (ignore isoform specification) and grab the mapping type and mapped value
        mapping_id = re.match("^([^-]*)", mapping_info[0]).groups()[0]
        mapping_type = mapping_info[1].lower()
        mapped_value = mapping_info[2]

        # If this ID is associated with both mapping types, then process it
        if mapping_id in paired_ids:
            if mapping_id not in id_to_labels:
                id_to_labels[mapping_id] = {"uniref": set(), "ortholog": set()}
            if mapping_type == args.uniref_version:
                id_to_labels[mapping_id]["uniref"].add(mapped_value)
            elif mapping_type == args.ortholog_type:
                id_to_labels[mapping_id]["ortholog"].add(mapped_value)

# Remove set of paired IDs to save memory
del paired_ids

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

# Open the output stream
# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# For each ortholog, write a line for each UniRef version ID that maps to it
for ortholog in ortholog_to_uniref.keys():
    for uniref in ortholog_to_uniref[ortholog]:
        output.write("\t".join([ortholog, uniref]) + os.linesep)

output.close()
