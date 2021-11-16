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

installation_directory = "{INSTALL_DIR}/"
"""
Directory where MetaLAFFA is installed. By default, this directory is used to store reference databases, supporting files, and as the original source of customizable configuration files when creating new MetaLAFFA projects.
"""

database_directory = installation_directory + "databases/"
"""
Directory containing databases to map reads to
"""

gene_normalization_directory = installation_directory + "gene_normalization_files/"
"""
Directory containing gene count normalization files
"""

gene_to_ortholog_directory = installation_directory + "gene_to_ortholog_maps/"
"""
Directory containing gene-to-ortholog mapping files
"""

ortholog_to_grouping_directory = installation_directory + "ortholog_to_grouping_maps/"
"""
Directory containing ortholog-to-grouping mapping files for aggregating ortholog abundances into functional groups
"""

log_directory = "logging/"
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

source_directory = installation_directory + "src/"
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

provenance_separator = "/"
"""
The character to use when separating steps in the pipeline in filenames, which are used to track provenance. If nested subdirectories are used, this will be the "/" character to change the single filename into nested subdirectories.
"""
if not nested_subdirectories:
    provenance_separator = "_"

required_reference_directories = [database_directory, gene_normalization_directory, gene_to_ortholog_directory, ortholog_to_grouping_directory, source_directory]
"""
List of directories for reference data that must exist for proper pipeline operation, these directories should be created during database preparation if they do not already exist.
"""

required_project_directories = [initial_data_directory, output_directory, benchmark_dir, log_directory, summary_directory, final_output_directory]
"""
List of directories for project-specific input, output, and logging data that must exist for proper pipeline operation, these directories should be created when creating a new project directory if they do not already exist.
"""