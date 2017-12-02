#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
#
# Author: Alex Eng
# Date: 12/1/2017

import numpy,pandas,argparse
from future import *

parser = argparse.ArgumentParser(description="Combines files containing rows of a summary table into a single table file")
parser.add_argument("summary_table_row_file", nargs="+", help="A row of the summary table (with column headers)")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Read the summary table rows
summary_table_rows = [pandas.read_table(x, sep="\t", header=0) for x in args.summary_table_row_file]

# Concatenate the rows
output_table = pandas.concat(summary_table_rows)

# Print the output
output_string = output_table.to_csv(args.output, sep="\t", header=True, index=False)

# If no output file was specified, then we print the table to standard output
if args.output == None:
	print(output_string)