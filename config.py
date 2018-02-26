# Input and output directories. All other directories are contained within these
output_dir="/net/borenstein/vol2/annotation_pipeline_testing/output/"
fastq_directory="/net/borenstein/vol2/annotation_pipeline_testing/data_test/"

#If you want the output split into nested subdiretories based on the different paramters. If False, then there is a single output directory with all the parameters concatenated together
many_subdirectories=False
#If you want to specific the samples to do in a text file. If None the code will determine all the samples that are present in the data directory
samples_oi=None
#Do you want to delete the intermediate files to save on disk space?
delete_intermediates=False

#Diamond Executable
aligner="/net/borenstein/vol1/PROGRAMS/diamond"
#Diamond computation options
memory=220
cpus=24
block_size=36
index_chunks=1
#Diamond running options
alignment_method="blastx"
sensitivity="" #Empty string for default (fast),  --sensitive or --more-sensitive to enable those options
top_percentage=1
max_e_value=0.001
db="/net/borenstein/vol1/DATA_DIAMONDDBs/KEGG/KEGG_7_15_2013/KEGG_gene_peptides.dmnd"

#Gene mapping options
count_method_gene="fractional" #fractional or whole
count_method_ko="fractional" #fractional or whole

#Hit filtering options
best_n_hits=10
filtering_method="best_hit" # best_hit, best_ko, best_N_hits, or best_N_kos

#Normalization options
norm_method="musicc"
musicc_correction_method="use_generic" # use_generic or learn_model

#ko functional summary options
mapping_matrix="/net/borenstein/vol1/DATA_REFERENCE/KEGG/KEGG_2013_07_15/KEGG_PARSED_2013_07_15/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15"
summary_method="fractional" # fractional or whole
functional_level="module" # module, pathway, bacterial_module, or bacterial_pathway

#Output directories
host_filtered_directory=output_dir + "1_host_filtered/"
duplicate_filtered_directory=output_dir + "2_duplicate_filtered/"
quality_filtered_directory=output_dir + "3_quality_filtered/"
diamond_output_directory=output_dir + "4_blast_results/"
diamond_filtered_directory=output_dir + "5_filtered_blast_results/"
diamond_counts_directory=output_dir + "6_gene_profiles/"
ko_counts_directory=output_dir + "7_ko_profiles/"
ko_normalized_directory=output_dir + "8_normalized_ko_profiles/"
module_profiles_directory=output_dir + "9_module_profiles/"
log_directory=output_dir + "logging/"
summary_directory=output_dir + "summaries/"
