#!/usr/bin/env python

import argparse
import sys
import os
import re
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Calculates ortholog abundances based on gene abundances")
parser.add_argument("gene_abundance_file", help="Gene abundance file to calculate ortholog abundances from")
parser.add_argument("counting_method", choices=["fractional", "whole"], help="The counting method to use")
parser.add_argument("gene_to_ortholog", help="File mapping genes to orthologs")
parser.add_argument("--ortholog_name", "-n", default="ortholog", help="The name of the ortholog system used (e.g. KO, COG, etc.), will be used as the column header for ortholog IDs (default: %(default)s)")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
	output = open(args.output, "w")

# Read associations between genes and orthologs
gene_to_ortholog_mapping = {}

# Add each ortholog that appears in the file to the set of orthologs associated with the indicated gene
f = custom_read(args.gene_to_ortholog)
for line in f:

	# Skip commented-out lines
	if line[0] != "#":

		# Ortholog IDs are the first column, gene IDs are the second column
		split_line = line.strip().split()
		ortholog = split_line[0]
		gene = split_line[1]

		# Strip prefix from the ortholog if present
		ortholog = re.sub("^.*:", "", ortholog)

		# If the gene isn't in the dictionary yet, give it a set of associated orthologs
		if gene not in gene_to_ortholog_mapping:
			gene_to_ortholog_mapping[gene] = set()

		# Add the gene-to-ortholog association to the dictionary
		gene_to_ortholog_mapping[gene].add(ortholog)

f.close()

# Initialize the variables used for keeping track of information while looping through each gene to calculate ortholog abundances
ortholog_abundances = {}

# Open the gene abundance file
f = custom_read(args.gene_abundance_file)

# Grab the sample name
sample_name = f.readline().strip().split()[1]

# Calculate ortholog abundances
for line in f:

	# Split the line into the gene and it's abundance
	gene_abundance_fields = line.strip().split()

	# Grab the gene and abundance
	gene = gene_abundance_fields[0]
	abundance = float(gene_abundance_fields[1])

	# If the gene is associated with any orthologs, calculate the ortholog abundance contributions for that gene
	if gene in gene_to_ortholog_mapping:

		contribution_factor = abundance

		# If we are using the fractional counting method, then we divide the abundance evenly among associated orthologs
		if args.counting_method == "fractional":
			contribution_factor /= float(len(gene_to_ortholog_mapping[gene]))

		# Add the contribution factor to the abundance of each associated ortholog
		for ortholog in gene_to_ortholog_mapping[gene]:

			# If the ortholog hasn't been seen before, initialize it's abundance at zero and then add to it
			if ortholog not in ortholog_abundances:
				ortholog_abundances[ortholog] = 0
			ortholog_abundances[ortholog] += contribution_factor

# Write the output table (two columns, ortholog name and ortholog abundance in the sample)
output.write("\t".join([args.ortholog_name, sample_name]) + os.linesep)
for ortholog in ortholog_abundances:
	output.write("\t".join([ortholog, str(ortholog_abundances[ortholog])]) + os.linesep)

output.close()
