#!/usr/bin/env python

import argparse
import sys
import os
from file_handling_lib import *

GENE_COLUMN_HEADER = "gene"

# Parse command line arguments
parser = argparse.ArgumentParser(description="Calculates gene abundances based on filtered blast hits")
parser.add_argument("filtered_blast_output_file", help="Filtered blast output file to calculate gene abundances from")
parser.add_argument("sample_name", help="Name of sample to use as column header for output gene abundances")
parser.add_argument("counting_method", choices=["fractional", "whole"], help="The counting method to use")
parser.add_argument("normalization_file", help="Path to file containing normalization factors for gene counts")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Parse gene normalization factors
gene_normalization_factors = {}
f = custom_read(args.normalization_file)
for line in f:

    # Skip commented-out lines
    if line[0] != "#":

        # Add each gene (first column) mapped to its normalization factor (second column)
        split_line = line.strip().split()
        gene_normalization_factors[split_line[0]] = float(split_line[2])

f.close()

# Initialize the variables used for keeping track of information while looping through each blast hit to calculate gene abundances
gene_abundances = {}
curr_read_name = None
curr_read_hits = []
# Calculate gene abundances (we assume all hits for the same read will be sequential)
f = custom_read(args.filtered_blast_output_file)
for line in f:

    # Split the line into the filtered blast output fields (read, hit, e value)
    blast_output_fields = line.strip().split()

    # Grab the read name and gene hit
    read_name = blast_output_fields[0]
    gene_hit = blast_output_fields[1]

    # If we've found a new read, calculate the gene abundance contributions for that read and reset the tracking variables
    if read_name != curr_read_name:

        # If the previous read had hits, calculate the gene abundance contributions for that read
        if len(curr_read_hits) > 0:

            partial_count_factor = 1

            # If we are using the fractional counting method, then we divide the count of 1 evenly among the hits
            if args.counting_method == "fractional":
                partial_count_factor = 1.0/float(len(curr_read_hits))

            # Add the partial count factor to the abundance of each gene hit
            for hit in curr_read_hits:

                # If the gene has not been seen before, initialize it's abundance to 0 then add the partial count factor
                if hit not in gene_abundances:
                    gene_abundances[hit] = 0
                gene_abundances[hit] += partial_count_factor

        # Reset the tracking variables
        curr_read_name = read_name
        curr_read_hits = []

    curr_read_hits.append(gene_hit)

# Don't forget to calculate gene abundance contributions for the last read
if len(curr_read_hits) > 0:

    partial_count_factor = 1

    # If we are using the fractional counting method, then we divide the count of 1 evenly among the hits
    if args.counting_method == "fractional":
        partial_count_factor = 1.0/float(len(curr_read_hits))

    # Add the partial count factor to the abundance of each gene hit
    for hit in curr_read_hits:

        # If the gene has not been seen before, initialize it's abundance to 0 then add the partial count factor
        if hit not in gene_abundances:
            gene_abundances[hit] = 0
        gene_abundances[hit] += partial_count_factor

# Normalize abundances by gene normalization factors
for gene in gene_abundances:
    gene_abundances[gene] /= float(gene_normalization_factors[gene])

# Write the output table (two columns, gene name and gene abundance in the sample)
output.write("\t".join([GENE_COLUMN_HEADER, args.sample_name]) + os.linesep)
for gene in gene_abundances:
    output.write("\t".join([gene, str(gene_abundances[gene])]) + os.linesep)
output.close()
