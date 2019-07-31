#!/net/borenstein/vol1/PROGRAMS/python2/bin/python
#
# Author: Alex Eng
# Date: 12/1/2017

SAMPLE_COLUMN_HEADER = "sample"

# Check why the path is being overwritten, but for now we make it look in our lab-controlled python lib directory
import sys
sys.path = ['', '/net/borenstein/vol1/PROGRAMS/python2/lib/python27.zip', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/plat-linux2', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-tk', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-old', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/lib-dynload', '/net/borenstein/vol1/PROGRAMS/python2/lib/python2.7/site-packages']
import numpy,pandas,argparse,re
from future import *

# Functions to parse file names for sample and file type information
def get_sample(filename):
	'''Returns the name of the sample associated with the input filename'''
	return re.match("^([^\.]*)\..*$", filename).group(1)

def get_file_type(filename):
	'''Returns the file type (R1, R2, S) of the input filename'''
	return re.match("^[^\.]*\.([^\.]*)\..*$", filename).group(1)

parser = argparse.ArgumentParser(description="Combines summary tables for each step of the pipeline into a single master summary table")
parser.add_argument("input_summary", help="The table summarizing information about the input files for the pipeline")
parser.add_argument("--host_filtering_summary", help="The table summarizing the host filtering step")
parser.add_argument("--duplicate_filtering_summary", help="The table summarizing the duplicate filtering step")
parser.add_argument("--quality_filtering_summary", help="The table summarizing the quality filtering step")
parser.add_argument("mapping_summary", help="The table summarizing the mapping step")
parser.add_argument("blast_hit_filtering_summary", help="The table summarizing the blast hit filtering step")
parser.add_argument("gene_counting_summary", help="The table summarizing the gene counting step")
parser.add_argument("ko_counting_summary", help="The table summarizing the ko counting step")
parser.add_argument("functional_level_summarization_summaries", nargs="*", help="The tables summarizing any functional level summarization steps")
parser.add_argument("--output", "-o", default=None, help="File to write output to (default: print to standard output)")
args = parser.parse_args()

####################################### Initialize master summary table #######################################
# Read the input summary file
input_summary = pandas.read_table(args.input_summary, sep="\t", header = 0)

# Add a column indicating the sample name for each file
input_summary[SAMPLE_COLUMN_HEADER] = input_summary["fastq_file"].apply(get_sample)

# Add a column indicating which file type (R1, R2, Singleton) each file is
input_summary["file_type"] = input_summary["fastq_file"].apply(get_file_type)

# Grab the subset file data for each file type for each sample
r1_input_summary = input_summary[input_summary["file_type"] == "R1"]
r2_input_summary = input_summary[input_summary["file_type"] == "R2"]
singleton_input_summary = input_summary[input_summary["file_type"] == "S"]

# Set the index for each subset table to be the sample
r1_input_summary = r1_input_summary.set_index(SAMPLE_COLUMN_HEADER)
r2_input_summary = r2_input_summary.set_index(SAMPLE_COLUMN_HEADER)
singleton_input_summary = singleton_input_summary.set_index(SAMPLE_COLUMN_HEADER)

# Rename file data columns to indicate the source file type
r1_input_summary = r1_input_summary.rename(index=str, columns={"reads":"r1_reads", "average_read_length":"r1_average_read_length", "average_base_quality":"r1_average_base_quality"})
r2_input_summary = r2_input_summary.rename(index=str, columns={"reads":"r2_reads", "average_read_length":"r2_average_read_length", "average_base_quality":"r2_average_base_quality"})
singleton_input_summary = singleton_input_summary.rename(index=str, columns={"reads":"singleton_reads", "average_read_length":"singleton_average_read_length", "average_base_quality":"singleton_average_base_quality"})

# Remove file name and file type columns
r1_input_summary = r1_input_summary.drop(columns=["fastq_file", "file_type"])
r2_input_summary = r2_input_summary.drop(columns=["fastq_file", "file_type"])
singleton_input_summary = singleton_input_summary.drop(columns=["fastq_file", "file_type"])

# Merge the input file summary tables into the initial master table
summary_table = pandas.concat([r1_input_summary, r2_input_summary, singleton_input_summary], axis=1)

# Calculate number of reads, average read length, and average base quality across all files for a single sample
summary_table_for_calculation = summary_table.fillna(0)

# Total reads is just the sum of reads from each file
summary_table["total_reads"] = summary_table_for_calculation["r1_reads"] + summary_table_for_calculation["r2_reads"] + summary_table_for_calculation["singleton_reads"]

# We first get the total number of bases per file (reads * average_read_length), add them together, then divide by the total number of reads
summary_table["overall_average_read_length"] = ((summary_table_for_calculation["r1_reads"] * summary_table_for_calculation["r1_average_read_length"]) + (summary_table_for_calculation["r2_reads"] * summary_table_for_calculation["r2_average_read_length"]) + (summary_table_for_calculation["singleton_reads"] * summary_table_for_calculation["singleton_average_read_length"]))/summary_table["total_reads"]

# We first get the total quality score of all bases per file (bases * average_base_quality), add them together, then divide by the number of bases
summary_table["overall_average_base_quality"] = ((summary_table_for_calculation["r1_reads"] * summary_table_for_calculation["r1_average_read_length"] * summary_table_for_calculation["r1_average_base_quality"]) + (summary_table_for_calculation["r2_reads"] * summary_table_for_calculation["r2_average_read_length"] * summary_table_for_calculation["r2_average_base_quality"]) + (summary_table_for_calculation["singleton_reads"] * summary_table_for_calculation["singleton_average_read_length"] * summary_table_for_calculation["singleton_average_base_quality"]))/(summary_table["total_reads"] * summary_table["overall_average_read_length"])

##################################### Incorporate host filtering summary #####################################
# If we performed host filtering and their is a host filtering summary, add that data to the summary table
if args.host_filtering_summary:

	# Read the host filtering summary file
	host_filtering_summary = pandas.read_table(args.host_filtering_summary, sep="\t", header=0)

	# Add a column indicating the sample name for each file
	host_filtering_summary[SAMPLE_COLUMN_HEADER] = host_filtering_summary["filtered_fastq_file"].apply(get_sample)

	# Add a column indicating which file type (R1, R2, Singleton) each file is
	host_filtering_summary["file_type"] = host_filtering_summary["filtered_fastq_file"].apply(get_file_type)

	# Grab the subset summary data for each file type for each sample
	r1_host_filtering_summary = host_filtering_summary[host_filtering_summary["file_type"] == "R1"]
	r2_host_filtering_summary = host_filtering_summary[host_filtering_summary["file_type"] == "R2"]
	singleton_host_filtering_summary = host_filtering_summary[host_filtering_summary["file_type"] == "S"]

	# Set the index for each subset table to be the sample
	r1_host_filtering_summary = r1_host_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)
	r2_host_filtering_summary = r2_host_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)
	singleton_host_filtering_summary = singleton_host_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)

	# Rename file data columns to indicate the source file type
	r1_host_filtering_summary = r1_host_filtering_summary.rename(index=str, columns={"post_host_filtering_reads":"r1_post_host_filtering_reads", "post_host_filtering_average_read_length":"r1_post_host_filtering_average_read_length", "post_host_filtering_average_base_quality":"r1_post_host_filtering_average_base_quality"})
	r2_host_filtering_summary = r2_host_filtering_summary.rename(index=str, columns={"post_host_filtering_reads":"r2_post_host_filtering_reads", "post_host_filtering_average_read_length":"r2_post_host_filtering_average_read_length", "post_host_filtering_average_base_quality":"r2_post_host_filtering_average_base_quality"})
	singleton_host_filtering_summary = singleton_host_filtering_summary.rename(index=str, columns={"post_host_filtering_reads":"singleton_post_host_filtering_reads", "post_host_filtering_average_read_length":"singleton_post_host_filtering_average_read_length", "post_host_filtering_average_base_quality":"singleton_post_host_filtering_average_base_quality"})

	# Remove file name and file type columns
	r1_host_filtering_summary = r1_host_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])
	r2_host_filtering_summary = r2_host_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])
	singleton_host_filtering_summary = singleton_host_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])

	# Merge the host filtering summary tables into the master table
	summary_table = pandas.concat([summary_table, r1_host_filtering_summary, r2_host_filtering_summary, singleton_host_filtering_summary], axis=1)

	# Calculate the number of reads, average read length, and average base quality across all files for a single sample post host filtering
	summary_table_for_calculation = summary_table.fillna(0)

	# Total reads is just the sum of reads from each file post host filtering
	summary_table["total_post_host_filtering_reads"] = summary_table_for_calculation["r1_post_host_filtering_reads"] + summary_table_for_calculation["r2_post_host_filtering_reads"] + summary_table_for_calculation["singleton_post_host_filtering_reads"]

	# We first get the total number of bases per file (reads * average_read_length), add them together, then divide by the total number of reads
	summary_table["overall_post_host_filtering_average_read_length"] = ((summary_table_for_calculation["r1_post_host_filtering_reads"] * summary_table_for_calculation["r1_post_host_filtering_average_read_length"]) + (summary_table_for_calculation["r2_post_host_filtering_reads"] * summary_table_for_calculation["r2_post_host_filtering_average_read_length"]) + (summary_table_for_calculation["singleton_post_host_filtering_reads"] * summary_table_for_calculation["singleton_post_host_filtering_average_read_length"]))/summary_table["total_post_host_filtering_reads"]

	# We first get the total quality score of all bases per file (bases * average_base_quality), add them together, then divide by the number of bases
	summary_table["overall_post_host_filtering_average_base_quality"] = ((summary_table_for_calculation["r1_post_host_filtering_reads"] * summary_table_for_calculation["r1_post_host_filtering_average_read_length"] * summary_table_for_calculation["r1_post_host_filtering_average_base_quality"]) + (summary_table_for_calculation["r2_post_host_filtering_reads"] * summary_table_for_calculation["r2_post_host_filtering_average_read_length"] * summary_table_for_calculation["r2_post_host_filtering_average_base_quality"]) + (summary_table_for_calculation["singleton_post_host_filtering_reads"] * summary_table_for_calculation["singleton_post_host_filtering_average_read_length"] * summary_table_for_calculation["singleton_post_host_filtering_average_base_quality"]))/(summary_table["total_post_host_filtering_reads"] * summary_table["overall_post_host_filtering_average_read_length"])

################################### Incorporate duplicate filtering summary ###################################
# If we performed duplicate filtering and their is a duplicate filtering summary, add that data to the summary table
if args.duplicate_filtering_summary:

	# Read the duplicate filtering summary file
	duplicate_filtering_summary = pandas.read_table(args.duplicate_filtering_summary, sep="\t", header=0)

	# Add a column indicating the sample name for each file
	duplicate_filtering_summary[SAMPLE_COLUMN_HEADER] = duplicate_filtering_summary["filtered_fastq_file"].apply(get_sample)

	# Add a column indicating which file type (R1, R2, Singleton) each file is
	duplicate_filtering_summary["file_type"] = duplicate_filtering_summary["filtered_fastq_file"].apply(get_file_type)

	# Grab the subset summary data for each file type for each sample
	r1_duplicate_filtering_summary = duplicate_filtering_summary[duplicate_filtering_summary["file_type"] == "R1"]
	r2_duplicate_filtering_summary = duplicate_filtering_summary[duplicate_filtering_summary["file_type"] == "R2"]
	singleton_duplicate_filtering_summary = duplicate_filtering_summary[duplicate_filtering_summary["file_type"] == "S"]

	# Set the index for each subset table to be the sample
	r1_duplicate_filtering_summary = r1_duplicate_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)
	r2_duplicate_filtering_summary = r2_duplicate_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)
	singleton_duplicate_filtering_summary = singleton_duplicate_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)

	# Rename file data columns to indicate the source file type
	r1_duplicate_filtering_summary = r1_duplicate_filtering_summary.rename(index=str, columns={"post_duplicate_filtering_reads":"r1_post_duplicate_filtering_reads", "post_duplicate_filtering_average_read_length":"r1_post_duplicate_filtering_average_read_length", "post_duplicate_filtering_average_base_quality":"r1_post_duplicate_filtering_average_base_quality"})
	r2_duplicate_filtering_summary = r2_duplicate_filtering_summary.rename(index=str, columns={"post_duplicate_filtering_reads":"r2_post_duplicate_filtering_reads", "post_duplicate_filtering_average_read_length":"r2_post_duplicate_filtering_average_read_length", "post_duplicate_filtering_average_base_quality":"r2_post_duplicate_filtering_average_base_quality"})
	singleton_duplicate_filtering_summary = singleton_duplicate_filtering_summary.rename(index=str, columns={"post_duplicate_filtering_reads":"singleton_post_duplicate_filtering_reads", "post_duplicate_filtering_average_read_length":"singleton_post_duplicate_filtering_average_read_length", "post_duplicate_filtering_average_base_quality":"singleton_post_duplicate_filtering_average_base_quality"})

	# Remove file name and file type columns
	r1_duplicate_filtering_summary = r1_duplicate_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])
	r2_duplicate_filtering_summary = r2_duplicate_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])
	singleton_duplicate_filtering_summary = singleton_duplicate_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])

	# Merge the duplicate filtering summary tables into the master table
	summary_table = pandas.concat([summary_table, r1_duplicate_filtering_summary, r2_duplicate_filtering_summary, singleton_duplicate_filtering_summary], axis=1)

	# Calculate the number of reads, average read length, and average base quality across all files for a single sample post duplicate filtering
	summary_table_for_calculation = summary_table.fillna(0)

	# Total reads is just the sum of reads from each file post duplicate filtering
	summary_table["total_post_duplicate_filtering_reads"] = summary_table_for_calculation["r1_post_duplicate_filtering_reads"] + summary_table_for_calculation["r2_post_duplicate_filtering_reads"] + summary_table_for_calculation["singleton_post_duplicate_filtering_reads"]

	# We first get the total number of bases per file (reads * average_read_length), add them together, then divide by the total number of reads
	summary_table["overall_post_duplicate_filtering_average_read_length"] = ((summary_table_for_calculation["r1_post_duplicate_filtering_reads"] * summary_table_for_calculation["r1_post_duplicate_filtering_average_read_length"]) + (summary_table_for_calculation["r2_post_duplicate_filtering_reads"] * summary_table_for_calculation["r2_post_duplicate_filtering_average_read_length"]) + (summary_table_for_calculation["singleton_post_duplicate_filtering_reads"] * summary_table_for_calculation["singleton_post_duplicate_filtering_average_read_length"]))/summary_table["total_post_duplicate_filtering_reads"]

	# We first get the total quality score of all bases per file (bases * average_base_quality), add them together, then divide by the number of bases
	summary_table["overall_post_duplicate_filtering_average_base_quality"] = ((summary_table_for_calculation["r1_post_duplicate_filtering_reads"] * summary_table_for_calculation["r1_post_duplicate_filtering_average_read_length"] * summary_table_for_calculation["r1_post_duplicate_filtering_average_base_quality"]) + (summary_table_for_calculation["r2_post_duplicate_filtering_reads"] * summary_table_for_calculation["r2_post_duplicate_filtering_average_read_length"] * summary_table_for_calculation["r2_post_duplicate_filtering_average_base_quality"]) + (summary_table_for_calculation["singleton_post_duplicate_filtering_reads"] * summary_table_for_calculation["singleton_post_duplicate_filtering_average_read_length"] * summary_table_for_calculation["singleton_post_duplicate_filtering_average_base_quality"]))/(summary_table["total_post_duplicate_filtering_reads"] * summary_table["overall_post_duplicate_filtering_average_read_length"])

#################################### Incorporate quality filtering summary ####################################
# If we performed quality filtering and their is a quality filtering summary, add that data to the summary table
if args.quality_filtering_summary:

	# Read the quality filtering summary file
	quality_filtering_summary = pandas.read_table(args.quality_filtering_summary, sep="\t", header=0)

	# Add a column indicating the sample name for each file
	quality_filtering_summary[SAMPLE_COLUMN_HEADER] = quality_filtering_summary["filtered_fastq_file"].apply(get_sample)

	# Add a column indicating which file type (R1, R2, Singleton) each file is
	quality_filtering_summary["file_type"] = quality_filtering_summary["filtered_fastq_file"].apply(get_file_type)

	# Grab the subset summary data for each file type for each sample
	r1_quality_filtering_summary = quality_filtering_summary[quality_filtering_summary["file_type"] == "R1"]
	r2_quality_filtering_summary = quality_filtering_summary[quality_filtering_summary["file_type"] == "R2"]
	singleton_quality_filtering_summary = quality_filtering_summary[quality_filtering_summary["file_type"] == "S"]

	# Set the index for each subset table to be the sample
	r1_quality_filtering_summary = r1_quality_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)
	r2_quality_filtering_summary = r2_quality_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)
	singleton_quality_filtering_summary = singleton_quality_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)

	# Rename file data columns to indicate the source file type
	r1_quality_filtering_summary = r1_quality_filtering_summary.rename(index=str, columns={"post_quality_filtering_paired_reads":"r1_post_quality_filtering_paired_reads", "post_quality_filtering_paired_average_read_length":"r1_post_quality_filtering_paired_average_read_length", "post_quality_filtering_paired_average_base_quality":"r1_post_quality_filtering_paired_average_base_quality", "post_quality_filtering_singleton_reads":"r1_post_quality_filtering_singleton_reads", "post_quality_filtering_singleton_average_read_length":"r1_post_quality_filtering_singleton_average_read_length", "post_quality_filtering_singleton_average_base_quality":"r1_post_quality_filtering_singleton_average_base_quality"})
	r2_quality_filtering_summary = r2_quality_filtering_summary.rename(index=str, columns={"post_quality_filtering_paired_reads":"r2_post_quality_filtering_paired_reads", "post_quality_filtering_paired_average_read_length":"r2_post_quality_filtering_paired_average_read_length", "post_quality_filtering_paired_average_base_quality":"r2_post_quality_filtering_paired_average_base_quality", "post_quality_filtering_singleton_reads":"r2_post_quality_filtering_singleton_reads", "post_quality_filtering_singleton_average_read_length":"r2_post_quality_filtering_singleton_average_read_length", "post_quality_filtering_singleton_average_base_quality":"r2_post_quality_filtering_singleton_average_base_quality"})
	singleton_quality_filtering_summary = singleton_quality_filtering_summary.rename(index=str, columns={"post_quality_filtering_paired_reads":"singleton_post_quality_filtering_paired_reads", "post_quality_filtering_paired_average_read_length":"singleton_post_quality_filtering_paired_average_read_length", "post_quality_filtering_paired_average_base_quality":"singleton_post_quality_filtering_paired_average_base_quality", "post_quality_filtering_singleton_reads":"singleton_post_quality_filtering_singleton_reads", "post_quality_filtering_singleton_average_read_length":"singleton_post_quality_filtering_singleton_average_read_length", "post_quality_filtering_singleton_average_base_quality":"singleton_post_quality_filtering_singleton_average_base_quality"})

	# Remove file name, file type, and singleton post quality filtering paired read columns
	r1_quality_filtering_summary = r1_quality_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])
	r2_quality_filtering_summary = r2_quality_filtering_summary.drop(columns=["filtered_fastq_file", "file_type"])
	singleton_quality_filtering_summary = singleton_quality_filtering_summary.drop(columns=["filtered_fastq_file", "file_type", "singleton_post_quality_filtering_paired_reads", "singleton_post_quality_filtering_paired_average_read_length", "singleton_post_quality_filtering_paired_average_base_quality"])

	# Merge the quality filtering summary tables into the master table
	summary_table = pandas.concat([summary_table, r1_quality_filtering_summary, r2_quality_filtering_summary, singleton_quality_filtering_summary], axis=1)

	# Calculate the number of reads, average read length, and average base quality for paired reads, singleton reads, and all reads across all files for a single sample post quality filtering
	summary_table_for_calculation = summary_table.fillna(0)

	# Total reads is just the sum of reads from each file post quality filtering
	summary_table["total_post_quality_filtering_paired_reads"] = summary_table_for_calculation["r1_post_quality_filtering_paired_reads"] + summary_table_for_calculation["r2_post_quality_filtering_paired_reads"]
	summary_table["total_post_quality_filtering_singleton_reads"] = summary_table_for_calculation["r1_post_quality_filtering_singleton_reads"] + summary_table_for_calculation["r2_post_quality_filtering_singleton_reads"] + summary_table_for_calculation["singleton_post_quality_filtering_singleton_reads"]
	summary_table["total_post_quality_filtering_reads"] = summary_table["total_post_quality_filtering_paired_reads"] + summary_table["total_post_quality_filtering_singleton_reads"]

	# We first get the total number of bases per file (reads * average_read_length), add them together, then divide by the total number of reads
	summary_table["overall_post_quality_filtering_paired_average_read_length"] = ((summary_table_for_calculation["r1_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r1_post_quality_filtering_paired_average_read_length"]) + (summary_table_for_calculation["r2_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r2_post_quality_filtering_paired_average_read_length"]))/summary_table["total_post_quality_filtering_paired_reads"]
	summary_table["overall_post_quality_filtering_singleton_average_read_length"] = ((summary_table_for_calculation["r1_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r1_post_quality_filtering_singleton_average_read_length"]) + (summary_table_for_calculation["r2_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r2_post_quality_filtering_singleton_average_read_length"]) + summary_table_for_calculation["singleton_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["singleton_post_quality_filtering_singleton_average_read_length"])/summary_table["total_post_quality_filtering_singleton_reads"]
	summary_table["overall_post_quality_filtering_average_read_length"] = ((summary_table_for_calculation["r1_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r1_post_quality_filtering_paired_average_read_length"]) + (summary_table_for_calculation["r1_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r1_post_quality_filtering_singleton_average_read_length"]) + (summary_table_for_calculation["r2_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r2_post_quality_filtering_paired_average_read_length"]) + (summary_table_for_calculation["r2_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r2_post_quality_filtering_singleton_average_read_length"]) + (summary_table_for_calculation["singleton_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["singleton_post_quality_filtering_singleton_average_read_length"]))/summary_table["total_post_quality_filtering_reads"]

	# We first get the total quality score of all bases per file (bases * average_base_quality), add them together, then divide by the number of bases
	summary_table["overall_post_quality_filtering_paired_average_base_quality"] = ((summary_table_for_calculation["r1_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r1_post_quality_filtering_paired_average_read_length"] * summary_table_for_calculation["r1_post_quality_filtering_paired_average_base_quality"]) + (summary_table_for_calculation["r2_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r2_post_quality_filtering_paired_average_read_length"] * summary_table_for_calculation["r2_post_quality_filtering_paired_average_base_quality"]))/(summary_table["total_post_quality_filtering_paired_reads"] * summary_table["overall_post_quality_filtering_paired_average_read_length"])
	summary_table["overall_post_quality_filtering_singleton_average_base_quality"] = ((summary_table_for_calculation["r1_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r1_post_quality_filtering_singleton_average_read_length"] * summary_table_for_calculation["r1_post_quality_filtering_singleton_average_base_quality"]) + (summary_table_for_calculation["r2_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r2_post_quality_filtering_singleton_average_read_length"] * summary_table_for_calculation["r2_post_quality_filtering_singleton_average_base_quality"]) + (summary_table_for_calculation["singleton_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["singleton_post_quality_filtering_singleton_average_read_length"] * summary_table_for_calculation["singleton_post_quality_filtering_singleton_average_base_quality"]))/(summary_table["total_post_quality_filtering_singleton_reads"] * summary_table["overall_post_quality_filtering_singleton_average_read_length"])
	summary_table["overall_post_quality_filtering_average_base_quality"] = ((summary_table_for_calculation["r1_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r1_post_quality_filtering_paired_average_read_length"] * summary_table_for_calculation["r1_post_quality_filtering_paired_average_base_quality"]) + (summary_table_for_calculation["r1_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r1_post_quality_filtering_singleton_average_read_length"] * summary_table_for_calculation["r1_post_quality_filtering_singleton_average_base_quality"]) + (summary_table_for_calculation["r2_post_quality_filtering_paired_reads"] * summary_table_for_calculation["r2_post_quality_filtering_paired_average_read_length"] * summary_table_for_calculation["r2_post_quality_filtering_paired_average_base_quality"]) + (summary_table_for_calculation["r2_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["r2_post_quality_filtering_singleton_average_read_length"] * summary_table_for_calculation["r2_post_quality_filtering_singleton_average_base_quality"]) + (summary_table_for_calculation["singleton_post_quality_filtering_singleton_reads"] * summary_table_for_calculation["singleton_post_quality_filtering_singleton_average_read_length"] * summary_table_for_calculation["singleton_post_quality_filtering_singleton_average_base_quality"]))/(summary_table["total_post_quality_filtering_reads"] * summary_table["overall_post_quality_filtering_average_read_length"])

######################################### Incorporate mapping summary #########################################
# Read the mapping summary table
mapping_summary = pandas.read_table(args.mapping_summary, sep="\t", header=0)

# Add a column indicating the sample name for each file
mapping_summary[SAMPLE_COLUMN_HEADER] = mapping_summary["blast_output_file"].apply(get_sample)

# Set the index to be the sample
mapping_summary = mapping_summary.set_index(SAMPLE_COLUMN_HEADER)

# Remove file name column
mapping_summary = mapping_summary.drop(columns=["blast_output_file"])

# Merge the mapping summary tables into the master table
summary_table = pandas.concat([summary_table, mapping_summary], axis=1)

################################### Incorporate blast hit filtering summary ###################################
# Read the blast hit filtering summary table
blast_hit_filtering_summary = pandas.read_table(args.blast_hit_filtering_summary, sep="\t", header=0)

# Add a column indicating the sample name for each file
blast_hit_filtering_summary[SAMPLE_COLUMN_HEADER] = blast_hit_filtering_summary["blast_output_file"].apply(get_sample)

# Set the index to be the sample
blast_hit_filtering_summary = blast_hit_filtering_summary.set_index(SAMPLE_COLUMN_HEADER)

# Rename file data columns to indicate that they are post hit filtering
blast_hit_filtering_summary = blast_hit_filtering_summary.rename(index=str, columns={"matched_reads":"post_hit_filtering_matched_reads", "matches":"post_hit_filtering_matches", "average_e_value":"post_hit_filtering_average_e_value"})

# Remove the file name column
blast_hit_filtering_summary = blast_hit_filtering_summary.drop(columns=["blast_output_file"])

# Merge the blast hit filtering summary table into the master table
summary_table = pandas.concat([summary_table, blast_hit_filtering_summary], axis=1)

###################################### Incorporate gene counting summary ######################################
# Read the gene counting summary table
gene_counting_summary = pandas.read_table(args.gene_counting_summary, sep="\t", header=0)

# Add a column indicating the sample name for each file
gene_counting_summary[SAMPLE_COLUMN_HEADER] = gene_counting_summary["gene_profile_file"].apply(get_sample)

# Set the index to be the sample
gene_counting_summary = gene_counting_summary.set_index(SAMPLE_COLUMN_HEADER)

# Remove the file name column
gene_counting_summary = gene_counting_summary.drop(columns=["gene_profile_file"])

# Merge the gene counting summary table into the master table
summary_table = pandas.concat([summary_table, gene_counting_summary], axis=1)

####################################### Incorporate ko counting summary #######################################
# Read the ko counting summary table
ko_counting_summary = pandas.read_table(args.ko_counting_summary, sep="\t", header=0)

# Add a column indicating the sample name for each file
ko_counting_summary[SAMPLE_COLUMN_HEADER] = ko_counting_summary["KO_profile_file"].apply(get_sample)

# Set the index to be the sample
ko_counting_summary = ko_counting_summary.set_index(SAMPLE_COLUMN_HEADER)

# Remove the file name column
ko_counting_summary = ko_counting_summary.drop(columns=["KO_profile_file"])

# Merge the ko counting summary table into the master table
summary_table = pandas.concat([summary_table, ko_counting_summary], axis=1)

############################ Incorporate functional level summarization summaries ############################
# Process the summary table for each functional level summarization performed
for filename in args.functional_level_summarization_summaries:

	# Read the functional level summarization summary tables
	functional_level_summarization_summary = pandas.read_table(filename, sep="\t", header=0)

	# Add a column indicating the sample name for each file
	functional_level_summarization_summary[SAMPLE_COLUMN_HEADER] = functional_level_summarization_summary["sample_file"]

	# Set the index to be the sample
	functional_level_summarization_summary = functional_level_summarization_summary.set_index(SAMPLE_COLUMN_HEADER)

	# Remove the file name column
	functional_level_summarization_summary = functional_level_summarization_summary.drop(columns=["sample_file"])

	# Merge the functional level summarization summary table into the master table
	summary_table = pandas.concat([summary_table, functional_level_summarization_summary], axis=1)

############################################## Write the output ##############################################
# Print the output
output_string = summary_table.to_csv(args.output, sep="\t", header=True, index=True, index_label = SAMPLE_COLUMN_HEADER)

# If no output file was specified, then we print the table to standard output
if args.output == None:
	print(output_string)