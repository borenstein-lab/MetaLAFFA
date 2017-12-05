#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
#
# Author: Alex Eng
# Date: 11/28/2017

KO_COLUMN_HEADER = "KO"

import argparse,sys
from file_handling import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Calculates KO abundances based on gene abundances")
parser.add_argument("gene_abundance_file", help="Gene abundance file to calculate KO abundances from")
parser.add_argument("counting_method", choices=["fractional", "whole"], help="The counting method to use")
parser.add_argument("--gene_to_ko_mapping", "-g", default="/net/borenstein/vol1/DATA_REFERENCE/KEGG/KEGG_2013_07_15/genes/ko/ko_genes.list", help="File mapping genes to KOs")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args=parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
	output = open(args.output, "w")

# Read associations between genes and kos
gene_to_ko_mapping = {}

# Add each ko that appears in the file to the set of kos associated with the indicated gene
f = custom_read(args.gene_to_ko_mapping)
for line in f:
	split_line = line.strip().split()
	ko = split_line[0]
	gene = split_line[1]

	# Strip the prefix from the ko if present
	if ko[:3] == "ko:":
		ko = ko[3:]

	# If the gene isn't in the dictionary yet, give it a set of associated kos
	if gene not in gene_to_ko_mapping:
		gene_to_ko_mapping[gene] = set()

	# Add the gene-to-ko association to the dicitonary
	gene_to_ko_mapping[gene].add(ko)

f.close()

# Initialize the variables used for keeping track of information while looping through each gene to calculate ko abundances
ko_abundances = {}

# Open the gene abundance file
f = custom_read(args.gene_abundance_file)

# Grab the sample name
sample_name = f.readline().strip().split()[1]

# Calculate ko abundances
for line in f:

	# Split the line into the gene and it's abundance
	gene_abundance_fields = line.strip().split()

	# Grab the gene and abundance
	gene = gene_abundance_fields[0]
	abundance = float(gene_abundance_fields[1])

	# If the gene is associated with any kos, calculate the ko abundance contributions for that gene
	if gene in gene_to_ko_mapping:

		contribution_factor = abundance

		# If we are using the fractional counting method, then we divide the abundance evenly among associated kos
		if args.counting_method == "fractional":
			contribution_factor /= float(len(gene_to_ko_mapping[gene]))

		# Add the contribution factor to the abundance of each associated ko
		for ko in gene_to_ko_mapping[gene]:

			# If the ko hasn't been seen before, initialize it's abundance at zero and then add to it
			if ko not in ko_abundances:
				ko_abundances[ko] = 0
			ko_abundances[ko] += contribution_factor

# Write the output table (two columns, ko name and ko abundance in the sample)
output.write("\t".join([KO_COLUMN_HEADER, sample_name]) + "\n")
for ko in ko_abundances:
	output.write("\t".join([ko, str(ko_abundances[ko])]) + "\n")
output.close()
