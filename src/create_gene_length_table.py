#!/usr/bin/env python

import argparse
import os
import sys
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Calculates gene lengths and the average gene length of a database of gene amino acid sequences")
parser.add_argument("fasta", help="FASTA file containing gene sequences")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args=parser.parse_args()

# Keep track of gene lengths (average length is printed at the top of the file, so we can't just stream to the output)
gene_lengths = {}

# Open the fasta
f = custom_read(args.fasta)

# Iterate through the file, keeping track of the name and length of the current gene, the total number of genes, and the total length
curr_name = ""
curr_length = 0
total_genes = 0
total_length = 0
for line in f:

	# If the line is not empty, process it
	if len(line) > 0:

		# If the line is a new gene, try to record information about the previous gene
		if line[0] == ">":

			#  If we have information about the previous gene, save it
			if curr_name != "" and curr_length > 0:
				gene_lengths[curr_name] = curr_length
				
			# Reset our tracking variables
			curr_name = line.strip().split()[0][1:]
			curr_length = 0

			# Increment the gene counter
			total_genes += 1

		# Otherwise, add the length of this line of sequence to the length of the gene and the total length
		else:
			curr_length += len(line.strip())
			total_length += len(line.strip())
f.close()

# If we have information about the previous gene, save it
if curr_name != "" and curr_length > 0:
	gene_lengths[curr_name] = curr_length

# Set output stream
output = sys.stdout
if args.output:
	output = open(args.output, "w")

# Calculate and write the average length as a comment at the top of the output
avg_length = float(total_length) / float(total_genes)
output.write("# " + str(avg_length) + os.linesep)

# Write the output table (three columns, gene name, gene length, and normalized gene length)
for gene in gene_lengths:
	output.write("\t".join([gene, str(gene_lengths[gene]), str(float(gene_lengths[gene]) / float(avg_length))]) + os.linesep)
output.close()
