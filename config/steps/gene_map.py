"""
Gene map step parameters
---------------------------

This configuration submodule contains parameters related to the gene map pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{wildcards.sample}.mapping.txt"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "gene_map"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.genes.tab"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {
    "memory": "40G"  # Custom memory request for gene mapping
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

required_programs = {}
"""
Dictionary defining the paths to programs used by this pipeline step
"""

non_essential_params = {}
"""
Dictionary defining the pipeline step's parameters that don't affect the output
"""

operating_params = {
    "type": "default",  # ID for operation to perform
    "method": "fractional"  # Gene counting method, defines how read counts are distributed among genes they map to, options include: fractional, whole
}
"""
Dictionary defining the pipeline step's parameters using when running the associated software
"""

benchmark_file = "{sample}.log"
"""
The benchmark filename pattern
"""

log_file = "{sample}.log"
"""
The log filename pattern
"""

# Defining options for different operations to run during this step


def default(inputs, outputs, wildcards, log):
    """
    Default gene map operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :param log: The log file
    :return: None.
    """

    # Set locations of reference files if necessary
    gene_normalization = op.gene_normalization_file

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        subprocess.run([fo.source_directory + "gene_map.py", inputs.input, wildcards.sample, operating_params["method"], gene_normalization, "--output", outputs[0]],
                       stdout=open(log[0], "a"),
                       stderr=subprocess.STDOUT)

    # Otherwise, if the input file is a dummy file, create a dummy output
    else:
        subprocess.run(["touch", outputs[0]])


# Defining the wrapper function that chooses which defined operation to run

def rule_function(inputs, outputs, wildcards, log):
    """
    How to run the software associated with this step

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :param log: The log file
    :return: None.
    """

    if operating_params["type"] == "default":
        default(inputs, outputs, wildcards, log)
