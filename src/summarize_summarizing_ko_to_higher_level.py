#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
#
# Author: Alex Eng
# Date: 11/30/2017

SAMPLE_NAME_COLUMN_HEADER = "sample_file"
FUNCTION_COUNT_COLUMN_HEADER = "number_of_functions"

# Check why the path is being overwritten, but for now we make it look in our lab-controlled python lib directory
import sys
sys.path = ['', '/net/borenstein/vol1/PROGRAMS/python2/lib/python27.zip', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/plat-linux2', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-tk', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-old', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-dynload', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/site-packages']
import numpy,pandas,argparse,os.path
from future import *

parser = argparse.ArgumentParser(description="Summarizes the results of summarizing KO profiles to a higher functional level")
parser.add_argument("profiles", help="The table of higher functional level profiles to summarize")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Read the profiles
profiles = pandas.read_table(args.profiles, sep = "\t", header = 0, index_col = 0)

# Count the number of non-zeros in each column (number of the higher level functions found in each sample)
output_table = profiles.astype(bool).sum(axis=0).to_frame()

# Rename the samples to be the base file name
output_table.index = [os.path.basename(x) for x in list(output_table.index)]

# Set name of function count column
output_table.columns = [FUNCTION_COUNT_COLUMN_HEADER]

# Print the output
output_string = output_table.to_csv(args.output, sep="\t", header=True, index=True, index_label = SAMPLE_NAME_COLUMN_HEADER)

# If no output file was specified, then we print the table to standard output
if args.output == None:
	print(output_string)
