"""
Map reads to genes step parameters
---------------------------

This configuration submodule contains parameters related to the map reads to genes pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{wildcards.sample}.{wildcards.type}.fastq"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "map_reads_to_genes"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.{type}.mapping.txt"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {
    "memory": "220G",  # Amount of RAM required for each read mapping cluster job
    "cores": 22,  # Number of cores requested for each read mapping cluster job
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

resource_params = {
    "diamond": "src/diamond-0.9.22/diamond",  # Location of the DIAMOND program
    "block_size": 36,  # Block size to use in DIAMOND
    "index_chunks": 1  # Index chunks for DIAMOND
}
"""
Dictionary defining the pipeline step's parameters that control resource usage but do not affect the output
"""

operating_params = {
    "type": "default",  # ID for operation to perform
    "method": "blastx",  # Alignment method
    "top_percentage": 1,  # Keep hits within X percent of the best score
    "evalue_cutoff": 0.001,  # Keep hits with e-value X or better
    "sensitivity": ""  # Sensitivity of DIAMOND mapping
}
"""
Dictionary defining the pipeline step's parameters using when running the associated software
"""

benchmark_file = "{sample}.{type}.log"
"""
The benchmark filename pattern
"""

# Defining options for different operations to run during this step


def default(inputs, outputs, wildcards):
    """
    Default FASTQ summary operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    # Set locations of reference files
    target_database = fo.database_directory + op.target_database + ".dmnd"

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        subprocess.run([resource_params["diamond"], operating_params["method"], "--block-size", str(resource_params["block_size"]), "--index-chunks", str(resource_params["index_chunks"]), "--threads", str(cluster_params["cores"] * 2), "--target_database", target_database, "--query", inputs.input, "--out", outputs[0], "--top", str(operating_params["top_percentage"]), "--evalue", str(operating_params["evalue_cutoff"]), operating_params["sensitivity"]])

    # Otherwise, if the input file is a dummy file, create a dummy output
    else:
        subprocess.run(["touch", outputs[0]])


# Defining the wrapper function that chooses which defined operation to run

def rule_function(inputs, outputs, wildcards):
    """
    How to run the software associated with this step

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    if operating_params["type"] == "default":
        default(inputs, outputs, wildcards)
