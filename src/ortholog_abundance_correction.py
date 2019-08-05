import pandas
import argparse

parser = argparse.ArgumentParser(description="Corrects ortholog abundances by normalizing sample ortholog counts to sum to 1")
parser.add_argument("ortholog_profiles", help="The table of ortholog profiles to aggregate")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Read the ortholog profiles and mapping matrix
ortholog_profiles = pandas.read_csv(args.ortholog_profiles, sep="\t", header=0, index_col=0)

# Replace NAs with zeros
ortholog_profiles = ortholog_profiles.fillna(0)

# Divide the ortholog counts in each sample by the sum of the ortholog counts for that sample
ortholog_profiles = ortholog_profiles.div(ortholog_profiles.sum(axis=0), axis=1)

# Print the output
output_string = ortholog_profiles.to_csv(args.output, sep="\t", header=True, index=True)

# If no output file was specified, then we print the table to standard output
if not args.output:
    print(output_string)
