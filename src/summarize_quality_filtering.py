#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
# Author: Alex Eng
# Date: 11/30/2017

NUM_LINES_PER_READ = 4

import argparse,sys
from file_handling import *
from future import *

# Parse command line arguments
parser = argparse.ArgumentParser(description="Summarizes quality read filtering")
parser.add_argument("original_r1_fastq", help="Pre-quality filtered R1 fastq file to determine the origin of new singletons")
parser.add_argument("original_r2_fastq", help="Pre-quality filtered R2 fastq file to determine the origin of new singletons")
parser.add_argument("filtered_r1_fastq", help="Filtered R1 fastq file for which to summarize quality filtering")
parser.add_argument("filtered_r2_fastq", help="Filtered R2 fastq file for which to summarize quality filtering")
parser.add_argument("filtered_new_singleton_fastq", help="Filtered fastq file of new singleton reads for which to summarize quality filtering")
parser.add_argument("filtered_old_singleton_fastq", help="Filtered fastq file of old singleton reads for which to summarize quality filtering")
parser.add_argument("--output", "-o", help="File to write output to (default: print to standard output)")
args = parser.parse_args()

# Set output stream
output = sys.stdout
if args.output:
    output = open(args.output, "w")

# Print the column names for the quality filtering summary file
output.write("\t".join(["filtered_fastq_file", "post_quality_filtering_paired_reads", "post_quality_filtering_singleton_reads"]) + "\n")

# Parse the new singleton file to get the read names of the new singletons
new_singletons = set()
f = custom_read(args.filtered_new_singleton_fastq)
line = f.readline()
count = 0
while line != "":

    # If we're on the first line of a read, then we found another read to add to the set of new singleton read names
    if count == 0:
        new_singletons.add(line.strip())

    # Now update our counter for the current line in the read's info
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()

# Parse the filtered R1 fastq to determine how many reads remain
r1_read_count = 0
f = custom_read(args.filtered_r1_fastq)
line = f.readline()
count = 0
while line != "":

    # If we're on the first line of a read, then we found another read
    if count == 0:
        r1_read_count += 1

    # Now update our counter for the current line in the read's info
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()

# Parse the filtered R2 fastq to determine how many reads remain
r2_read_count = 0
f = custom_read(args.filtered_r2_fastq)
line = f.readline()
count = 0
while line != "":

    # If we're on the first line of a read, then we found another read
    if count == 0:
        r2_read_count += 1

    # Now update our counter for the current line in the read's info
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()

# Parse the filtered singleton fastq to determine how many reads remain
singleton_read_count = 0
f = custom_read(args.filtered_old_singleton_fastq)
line = f.readline()
count = 0
while line != "":

    # If we're on the first line of a read, then we found another read
    if count == 0:
        singleton_read_count += 1

    # Now update our counter for the current line in the read's info
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()

# Parse the unfiltered R1 fastq to determine how many new singletons came from the R1 fastq
r1_new_singleton_read_count = 0
f = custom_read(args.original_r1_fastq)
line = f.readline()
count = 0
while line != "":

    # If we're on the first line of a read, then check it's name against the new singletons
    if count == 0:
        if line.strip() in new_singletons:
            r1_new_singleton_read_count += 1

    # Now update our counter for the current line in the read's info
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()

# Parse the unfiltered R2 fastq to determine how many new singletons came from the R2 fastq
r2_new_singleton_read_count = 0
f = custom_read(args.original_r2_fastq)
line = f.readline()
count = 0
while line != "":

    # If we're on the first line of a read, then check it's name against the new singletons
    if count == 0:
        if line.strip() in new_singletons:
            r2_new_singleton_read_count += 1

    # Now update our counter for the current line in the read's info
    count += 1
    count = count % NUM_LINES_PER_READ
    line = f.readline()
f.close()

# Print the quality filtering summary for this sample
output.write("\t".join([args.filtered_r1_fastq, str(r1_read_count), str(r1_new_singleton_read_count)]) + "\n")
output.write("\t".join([args.filtered_r2_fastq, str(r2_read_count), str(r2_new_singleton_read_count)]) + "\n")
output.write("\t".join([args.filtered_old_singleton_fastq, "0", str(singleton_read_count)]) + "\n")