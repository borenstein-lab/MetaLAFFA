#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
#
# Author: Alex Eng
# Date: 11/30/2017

FUNCTION_COLUMN_HEADER = "Function"

# Check why the path is being overwritten, but for now we make it look in our lab-controlled python lib directory
import sys
sys.path = ['', '/net/borenstein/vol1/PROGRAMS/python2/lib/python27.zip', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/plat-linux2', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-tk', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-old', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-dynload', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/site-packages']
import numpy,pandas,argparse
from future import *

parser = argparse.ArgumentParser(description="Summarizes KO profiles to a higher functional level")
parser.add_argument("ko_profiles", help="The table of KO profiles to summarize")
parser.add_argument("summary_method", choices=["fractional", "whole"], help="The method to use to map KOs to the higher functional level")
parser.add_argument("mapping_matrix", help="The matrix that maps KOs to the higher functional levels they belong to")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Read the KO profiles and mapping matrix
ko_profiles = pandas.read_table(args.ko_profiles, sep = "\t", header = 0, index_col = 0)
mapping_matrix = pandas.read_table(args.mapping_matrix, sep="\t", header=0, index_col = 0)

# Filter and sort the mapping matrix
mapping_matrix = mapping_matrix.reindex(list(ko_profiles.index))

# If using the fractional summary method, we divide KO contributions to higher functional levels evenly among the higher functions each KO belongs to
if args.summary_method == "fractional":
	
	# Divide each row by the sum of the row
	mapping_matrix = mapping_matrix.div(mapping_matrix.sum(axis=1), axis=0)
	mapping_matrix = mapping_matrix.fillna(0)

# Transpose the mapping matrix
mapping_matrix = mapping_matrix.transpose()	

# Multiply the mapping matrix and KO profiles
output_table = mapping_matrix.dot(ko_profiles)

# Remove any row that is all zeros
output_table = output_table[(output_table.T != 0).any()]

# Print the output
output_string = output_table.to_csv(args.output, sep="\t", header=True, index=True, index_label = FUNCTION_COLUMN_HEADER)

# If no output file was specified, then we print the table to standard output
if args.output == None:
	print(output_string)