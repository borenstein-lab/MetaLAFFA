#!/usr/bin/env python

import argparse
import sys
import os
from file_handling_lib import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Filters blast hits based on the selected criterion")
parser.add_argument("blast_output_file", help="blast output file to filter")
parser.add_argument("filtering_method", choices=["best_hit", "best_ortholog", "best_N_hits", "best_N_orthologs"], help="The filtering criterion to use")
parser.add_argument("--number", "-n", type=int, default=20, help="The number of best hits or kos to use")
parser.add_argument("--gene_to_ortholog_map", "-g", help="Path to file mapping genes to orthologs")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
	output = open(args.output, "w")

# If we are using an ortholog-based criterion, get the set of genes with associated orthologs
ortholog_genes = set()
if args.filtering_method in ["best_ortholog", "best_N_orthologs"]:

	# Check if a gene-to-ortholog mapping was provided, exit and report if not
	if not args.gene_to_ortholog_map:
		sys.exit("Filter method %s requires a gene-to-ortholog mapping file and none was provided." % args.filtering_method)

	# Add each gene that appears in the file to the set of genes with ortholog associations (second column)
	f = custom_read(args.args.gene_to_ortholog_map)
	for line in f:
		ortholog_genes.add(line.strip().split()[1])
	f.close()

# Initialize variables used for filtering blast hits
curr_read_name = None
curr_read_hits = []
curr_best_e_val = float("Inf")
curr_num_hits = 0

# Filter the blast output (we assume all hits for the same read will be sequential in ascending order of e value)
f = custom_read(args.blast_output_file)
for line in f:

	# Split the line into the blast output fields
	blast_output_fields = line.strip().split()

	# Grab the read name, gene hit, and e value
	read_name = blast_output_fields[0]
	gene_hit = blast_output_fields[1]
	e_val = float(blast_output_fields[10])

	# If we've found a new read, print any filtered hit information for the previous read and reset the tracking variables
	if read_name != curr_read_name:

		# Print the previous read hit information
		for hit in curr_read_hits:
			output.write(hit + os.linesep)

		# Reset the tracking variables for the new current read
		curr_read_name = read_name
		curr_read_hits = []
		curr_best_e_val = float("Inf")
		curr_num_hits = 0

	# If the current e value is less than or equal to the best e value for hits for this read and we're using a best hit method or we're using a best ortholog method and the gene is associated with a ortholog then we keep this hit
	if e_val <= curr_best_e_val and (args.filtering_method in ["best_hit", "best_N_hits"] or (args.filtering_method in ["best_ortholog", "best_N_orthologs"] and gene_hit in ortholog_genes)):
		curr_read_hits.append("\t".join([read_name, gene_hit, str(e_val)]))
		curr_num_hits += 1

		# If the e value was less than the best e value for this read, that means we had not found a read that satisfies the criterion requirements for this read yet and we should update the best e value for this read
		if e_val < curr_best_e_val:
			curr_best_e_val = e_val

	# Otherwise, if the current e value is greater than the best e value, then we only keep it if we are using a best N method, we haven't seen N hits yet, and if we are using the best N ortholog method then we also check if the gene is associated with aan ortholog
	elif e_val > curr_best_e_val and (args.filtering_method == "best_N_hits" or (args.filtering_method == "best_N_orthologs" and gene_hit in ortholog_genes)) and curr_num_hits < args.number:
		curr_read_hits.append("\t".join([read_name, gene_hit, str(e_val)]))
		curr_num_hits += 1

# Don't forget to print any info for the last read
for hit in curr_read_hits:
	output.write(hit + os.linesep)
output.close()
