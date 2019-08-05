"""
File organization parameters
----------------------------

This configuration submodule contains parameters related to organization of the pipeline file structure.
"""

initial_data_directory = "data/"
"""
Directory containing the input FASTQ files. The default assumes that you have created a local directory called *data/* where the FASTQs are located.
"""

output_directory = "output/"
"""
Directory containing intermediate and final output files. This directory will be created when the pipeline is run for the first time.
"""

benchmark_dir = "benchmark/"
"""
Directory containing benchmarking files. This directory will be created when the pipeline is run for the first time.
"""

bitmask_directory = "bitmasks/"
"""
Directory containing database bitmaps
"""

database_directory = "databases/"
"""
Directory containing databases to map reads to
"""

index_directory = "indices/"
"""
Directory containing indices for srprism
"""

gene_normalization_directory = "gene_normalization_files/"
"""
Directory containing gene count normalization files
"""

gene_to_ortholog_directory = "gene_to_ortholog_maps/"
"""
Directory containing gene-to-ortholog mapping files
"""

ortholog_to_grouping_directory = "ortholog_to_grouping_maps/"
"""
Directory containing ortholog-to-grouping mapping files for aggregating ortholog abundances into functional groups
"""

log_directory = output_directory + "logging/"
"""
Output subdirectory for pipeline log files.
"""

summary_directory = output_directory + "summaries/"
"""
Output subdirectory for pipeline step summary files.
"""

final_output_directory = output_directory + "final/"
"""
Output subdirectory in which to collect the final output files.
"""

source_directory = "src/"
"""
Directory for scripts and third-party tools used in the pipeline.
"""

tmp_dir = "/tmp/"
"""
Directory containing temporary intermediate files. These tend to be large uncompressed intermediate files that we need to compress before saving, if we are saving intermediate files.
"""

nested_subdirectories = True
"""
Whether provenance for the output for each step should be tracked in a nested subdirectory structure with subdirectory titles indicating parameters for previous processing steps, or whether previous processing step parameters should be encoded in the name of the output file. This can be used when the filesystem has limits on file name length.
"""

provenance_separator = "_"
"""
The character to use when separating steps in the pipeline in filenames, which are used to track provenance. If nested subdirectories are used, this is automatically changed to a "/" character to change the single filename into nested subdirectories.
"""
if nested_subdirectories:
    provenance_separator = "/"

required_directories = [initial_data_directory, output_directory, benchmark_dir, bitmask_directory, database_directory, index_directory, gene_normalization_directory, gene_to_ortholog_directory, ortholog_to_grouping_directory, log_directory, summary_directory, final_output_directory, source_directory]
"""
List of directories that must exist for proper pipeline operation, these directories should be created during pipeline setup if they do not already exist.
"""
