"""
File organization parameters
----------------------------

This configuration submodule contains parameters related to organization of the pipeline file structure.
"""

fastq_directory = "data/"
"""
Directory containing the input FASTQ files. The default assumes that you have created a local directory called *data/* where the FASTQs are located.
"""

output_dir = "output/"
"""
Directory containing intermediate and final output files. This directory will be created when the pipeline is run for the first time.
"""

host_filtered_fastq_directory = output_dir + "1_host_filtered/"
"""
Output subdirectory for host filtering FASTQs.
"""

duplicate_filtered_fastq_directory = output_dir + "2_duplicate_filtered/"
"""
Output subdirectory for duplicate filtered FASTQs.
"""

quality_filtered_fastq_directory = output_dir + "3_quality_filtered/"
"""
Output subdirectory for quality filtered FASTQs.
"""

alignment_directory = output_dir + "4_blast_results/"
"""
Output subdirectory for read alignments.
"""

filtered_alignment_directory = output_dir + "5_filtered_blast_results/"
"""
Output subdirectory for filtered read alignments.
"""

gene_profile_directory = output_dir + "6_gene_profiles/"
"""
Output subdirectory for gene profiles.
"""

ortholog_profile_directory = output_dir + "7_ko_profiles/"
"""
Output subdirectory for ortholog profiles.
"""

normalized_ortholog_profile_directory = output_dir + "8_normalized_ko_profiles/"
"""
Output subdirectory for normalized ortholog profiles.
"""

summarized_profile_directory = output_dir + "9_module_profiles/"
"""
Output subdirectory for functional profiles summarized to a higher functional classification level.
"""

log_directory = output_dir + "logging/"
"""
Output subdirectory for pipeline log files.
"""

summary_directory = output_dir + "summaries/"
"""
Output subdirectory for pipeline step summary files.
"""

tmp_dir = "/tmp/"
"""
Directory containing temporary intermediate files. These tend to be large uncompressed intermediate files that we need to compress before saving, if we are saving intermediate files.
"""

nested_subdirectories = False
"""
Whether provenance for the output for each step should be tracked in a nested subdirectory structure with subdirectory titles indicating parameters for previous processing steps, or whether previous processing step parameters should be encoded in the name of the output file. This can be used when the filesystem has limits on file name length.
"""