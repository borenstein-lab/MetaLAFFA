"""
Pipeline operation parameters
-----------------------------

This configuration submodule contains parameters related to user operation of the pipeline.
"""

import config.file_organization as fo

pipeline_step_list = "pipeline_steps.txt"
"""
The file listing which pipeline steps to run and defining how the inputs and outputs of different steps are linked. Each line of the file should consist of the name of a pipeline step (the name of a snakemake rule), followed by a colon and then an ordered comma-separated list of pipeline steps that are the input to the given step. Input steps should be ordered based on their usage in the step they feed into. For example:

step3:step1,step2

indicates that the pipeline step named "step3" should be run, and it requires the outputs of the "step1" and "step2" pipeline steps. 

There are five special cases in pipeline step definition:

The first is for steps that take initial input to the pipeline, which should have this listed as "INPUT" in the input steps (i.e. "step1:INPUT"). In this case, the pipeline will look for input files that match the step's default input file name format.

The second is for steps that produce desired final output files. These steps should be marked by including a "$" at the beginning of the step name.

The third is for steps that produce a summary what happened during a pipeline step. These steps should be marked by including a "*" at the beginning of the step name.

The fourth is for steps that take as input all pipeline step summaries (the output of steps marked with a "*"). The input to these steps should be defined as "SUMMARIES".

The fifth is for steps that should have their intermediate outputs deleted once they are no longer needed. These steps should be marked by including a "~" at the beginning of the step name.  

Steps can be skipped by commenting them out with the "#" character, though note that any later steps relying on skipped steps should have their input steps updated.
"""

conda_env = "{CONDA_ENV}"
"""
The name of the conda environment in which MetaLAFFA should be run.
"""

snakefile = "Snakefile"
"""
The file that defines what pipeline steps are available to Snakemake.
"""

sample_list = None
"""
The file to specifies which samples to annotate. If "None", then samples will be determined based on the contents of the input FASTQ directory.
"""

dummy_input_log = "dummy_input_files_generated.txt"
"""
The log file for saving which dummy files were generated because a sample had some, but not all, input files present.
"""

work_in_tmp_dir = True
"""
Whether outputs should first be generated in a temporary working directory before being moved to their final output directory. This is most helpful when you have local disk space constraints and need to zip output files because a raw output file can be generated in a cluster node's potentially larger temporary local directory, zipped, and then the zipped version of the output file can be moved to the final output location.
"""

zipped_files = True
"""
Whether input and output files are expected to be zipped. If true, intermediate files will always be zipped after creation.
"""

delete_intermediates = False
"""
Whether output files from intermediate pipeline steps should be saved. If not, then they are deleted after they are no longer required.
"""

fastq_types = {
    "forward": "R1",
    "reverse": "R2",
    "singleton": "S"
}
"""
Define the FASTQ type markers using in file naming.
"""

cpu_to_thread_multiplier = 2
"""
Define how many threads there are per cpu
"""

host_database = "hs37d5"
"""
The name of the host database to map reads to for host filtering
"""

target_database = "uniref90"
"""
The name of the target database to map reads to for gene identification 
"""

target_ortholog = "ko"
"""
The name of the target ortholog to map genes to for functional annotation
"""

gene_to_ortholog_mapping = target_database + "_to_" + target_ortholog
"""
The mapping to use for linking genes to orthologs
"""

ortholog_to_grouping_mappings = ["KOvsMODULE_BACTERIAL_KEGG_2013_07_15", "KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15"]
"""
List of ortholog-to-grouping mapping files for aggregating ortholog abundances into functional groups 
"""

wildcard_restrictions = {
    "type": list(fastq_types.values()),
    "mapping": ortholog_to_grouping_mappings
}
"""
Define restrictions on wildcard values for file name patterns
"""

snakemake = "snakemake"
"""
The path to the snakemake executable to use when running the pipeline
"""

bowtie2_build = "bowtie2-build"
"""
The path to the bowtie2-build executable to use when processing default databases
"""

diamond_db_suffix = ".dmnd"
"""
Suffix for a DIAMOND database
"""

gene_normalization_suffix = ".norm"
"""
Suffix for a gene count normalization file
"""

gene_to_ortholog_suffix = ".map"
"""
Suffix for a gene-to-ortholog mapping file
"""

ortholog_to_grouping_suffix = ".map"
"""
Suffix for an ortholog-to-grouping mapping file
"""

host_database_file = fo.database_directory + host_database + ".fa"
"""
The path to the database for host filtering
"""

host_index = fo.database_directory + host_database
"""
The path to the index file for host filtering
"""

target_database_file = fo.database_directory + target_database + diamond_db_suffix
"""
The path to the database for mapping reads to gene counts
"""

gene_normalization_file = fo.gene_normalization_directory + target_database + gene_normalization_suffix
"""
The path to the gene normalization file for normalization mapped gene counts
"""

gene_to_ortholog_file = fo.gene_to_ortholog_directory + gene_to_ortholog_mapping + gene_to_ortholog_suffix
"""
The path to the gene-to-ortholog mapping file for mapping gene counts to ortholog counts
"""

ortholog_to_grouping_files = [fo.ortholog_to_grouping_directory + mapping + ortholog_to_grouping_suffix for mapping in ortholog_to_grouping_mappings]
"""
The paths to the ortholog-to-grouping mapping files for mapping ortholog counts to higher functional grouping counts
"""
