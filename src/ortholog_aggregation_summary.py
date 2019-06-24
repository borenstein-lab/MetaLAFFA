import pandas
import argparse

SAMPLE_COLUMN_HEADER = "sample"

parser = argparse.ArgumentParser(description="Summarizes the results of summarizing ortholog profiles to a higher functional grouping")
parser.add_argument("profiles", help="The table of higher functional grouping profiles to summarize")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Read the profiles
profiles = pandas.read_csv(args.profiles, sep="\t", header=0, index_col=0)

# Get the name of the functional groupo that was summarized to
functional_group = profiles.index.name

# Count the number of non-zeros in each column (number of the higher groupings found in each sample)
output_table = profiles.astype(bool).sum(axis=0).to_frame()

# Set name of function count column
output_table.columns = [functional_group]

# Print the output
output_string = output_table.to_csv(args.output, sep="\t", header=True, index=True, index_label=SAMPLE_COLUMN_HEADER)

# If no output file was specified, then we print the table to standard output
if not args.output:
	print(output_string)
