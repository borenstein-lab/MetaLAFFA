#sample_file_name="/net/borenstein/vol2/annotation_pipeline_testing/samples.tab"
sample_file_name="samples.tab"
fastq_directory="/net/borenstein/vol2/annotation_pipeline_testing/data/"
diamond_output_directory="/net/borenstein/vol2/annotation_pipeline_testing/output/blast_results/"
log_directory="/net/borenstein/vol2/annotation_pipeline_testing/logging/"
output_directory="/net/borenstein/vol2/annotation_pipeline_testing/output/"
memory=220
cpus=24
block_size=36
index_chunks=1
aligner="/net/borenstein/vol1/PROGRAMS/diamond"
alignment_method="blastx"
sensitivity=None
top_percentage=1
max_e_value=0.001
count_method="fractional"
best_n_hits=10
filtering_method="best_hit"
bash_wrapper="/net/borenstein/vol1/PIPELINE/CODE_GENERAL/qsub_wrapper.sh"
#diamond_wrapper="/net/borenstein/vol2/AlexE/diamond_wrapper.sh"
#summarize_mapping="/net/borenstein/vol2/AlexE/summarize_mapping.py"
#summarize_mapping_wrapper="/net/borenstein/vol2/AlexE/summarize_mapping_wrapper.sh"
db="/net/borenstein/vol1/DATA_DIAMONDDBs/KEGG/KEGG_7_15_2013/KEGG_gene_peptides.dmnd"