import pandas
import argparse

parser = argparse.ArgumentParser(description="Aggregates ortholog profiles to a higher functional grouping")
parser.add_argument("ortholog_profiles", help="The table of ortholog profiles to aggregate")
parser.add_argument("summary_method", choices=["fractional", "whole"], help="The method to use to map orthologs to the higher functional grouping")
parser.add_argument("ortholog_to_grouping", help="The file mapping orthologs to groups defined at a higher functional level")
parser.add_argument("--grouping_name", "-n", default="group", help="The name of the grouping that orthologs are being aggregated into (e.g. module, pathway, etc.), will be used as the column header for grouping IDs (default: %(default)s)")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Read the ortholog profiles and mapping matrix
ortholog_profiles = pandas.read_csv(args.ortholog_profiles, sep="\t", header=0, index_col=0)
mapping_matrix = pandas.read_csv(args.ortholog_to_grouping, sep="\t", header=0, index_col=0)

# Replace NAs with zeros
ortholog_profiles = ortholog_profiles.fillna(0)

# Filter down to orthologs that are both in the profiles and in the mapping matrix
profile_orthologs = set(ortholog_profiles.index)
mapping_orthologs = set(mapping_matrix.index)
shared_orthologs = list(profile_orthologs.intersection(mapping_orthologs))

# Filter and sort the matrices
ortholog_profiles = ortholog_profiles.reindex(shared_orthologs)
mapping_matrix = mapping_matrix.reindex(shared_orthologs)

# If using the fractional summary method, we divide ortholog contributions to higher functional groups evenly among the higher groups each ortholog belongs to
if args.summary_method == "fractional":

    # Divide each row by the sum of the row
    mapping_matrix = mapping_matrix.div(mapping_matrix.sum(axis=1), axis=0)
    mapping_matrix = mapping_matrix.fillna(0)

# Transpose the mapping matrix
mapping_matrix = mapping_matrix.transpose()

# Multiply the mapping matrix and ortholog profiles
output_table = mapping_matrix.dot(ortholog_profiles)

# Remove any row that is all zeros
output_table = output_table[(output_table.T != 0).any()]

# Print the output
output_string = output_table.to_csv(args.output, sep="\t", header=True, index=True, index_label=args.grouping_name)

# If no output file was specified, then we print the table to standard output
if not args.output:
    print(output_string)
