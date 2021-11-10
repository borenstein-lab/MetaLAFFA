#!/usr/bin/env python

import argparse
import sys
import os
import re
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(
    description="Using a UniRef ID mapping file, creates a table mapping UniRef gene IDs to ortholog groups.")
parser.add_argument("uniref_to_uniprot", help="File from UniProt mapping UniRef IDs to UniProt IDs")
parser.add_argument("uniprot_to_ortholog", help="File from UniProt mapping UniProt IDs to ortholog IDs")
parser.add_argument("ortholog_type", help="The identifier for the type of ortholog to map genes to (e.g. ko, eggnog)")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Standardize strings to be lower-case
args.ortholog_type = args.ortholog_type.lower()

uniprot_to_ortholog = {}
with custom_read(args.uniprot_to_ortholog) as ortholog_mapping:
    curr_accession = None
    curr_orthologs = set()
    for line in ortholog_mapping:

        fields = line.strip().split()

        # If we hit a new accession ID, add any available accession data to the mapping dictionary and reset trackers
        # Trailing semicolons are stripped from all fields except first
        if fields[0] == "AC" and len(fields) > 1:
            if len(curr_orthologs) > 0:
                uniprot_to_ortholog[curr_accession] = curr_orthologs
            curr_accession = fields[1][:-1]
            curr_orthologs = set()

        # Otherwise, if we find a mapping to the ortholog type of interest, add that to the set of orthologs associated with our current accession ID
        elif fields[0] == "DR" and len(fields) > 2 and fields[1][:-1].lower() == args.ortholog_type:
            curr_orthologs.add(fields[2][:-1])

    # Don't forget last entry
    if len(curr_orthologs) > 0:
        uniprot_to_ortholog[curr_accession] = curr_orthologs

# Open the output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Scan UniRef-to-UniProt mapping file and write ortholog-to-UniRef mappings to output file
with custom_read(args.uniref_to_uniprot) as uniref_mapping:
    curr_uniref = None
    curr_orthologs = set()
    for line in uniref_mapping:

        fields = line.strip().split()

        # If we hit a new UniRef ID, write out ortholog-to-UniRef mappings for this UniRef ID
        if fields[0] != curr_uniref:
            if len(curr_orthologs) > 0:
                curr_orthologs = list(curr_orthologs)
                for ortholog in curr_orthologs:
                    output.write("\t".join([ortholog, curr_uniref]) + os.linesep)

            curr_uniref = fields[0]
            curr_orthologs = set()

        curr_uniprot = fields[1]
        if curr_uniprot in uniprot_to_ortholog:
            curr_orthologs = set.union(curr_orthologs, uniprot_to_ortholog[curr_uniprot])

output.close()
