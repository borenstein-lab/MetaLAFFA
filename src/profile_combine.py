import pandas
import pandas.errors
import argparse
import sys
import os

parser = argparse.ArgumentParser(description="Combines separate profile tables into a single profile table")
parser.add_argument("profile_tables", nargs="+", help="The profile tables to combine")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Read the summary tables
tables = []
for table in args.profile_tables:
    try:
        table = pandas.read_csv(table, sep="\t", header=0)
        tables.append(table)
    except pandas.errors.EmptyDataError:
        sys.stderr.write("Skipping empty table: " + table + os.linesep)

# Get the ID column name
id_column = tables[0].columns[0]

# Set indices of the tables
for table_index in range(len(tables)):
    tables[table_index] = tables[table_index].set_index(id_column)

# Group tables by the columns they share for concatenation
tables_by_columns = {}
for table in tables:
    table_cols = frozenset(table.columns)
    if table_cols not in tables_by_columns:
        tables_by_columns[table_cols] = []
    tables_by_columns[table_cols].append(table)

# Concatenate tables with shared column names
concatenated_tables = []
for table_cols in tables_by_columns:
    concatenated_tables.append(pandas.concat(tables_by_columns[table_cols]))

# Next concatenate tables by the index to get final table and write the table
final_table = pandas.concat(concatenated_tables, axis=1)
final_table.to_csv(path_or_buf=output, sep="\t", header=True, index=True, index_label=id_column, na_rep="0")

output.close()
