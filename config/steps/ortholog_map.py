"""
Ortholog map step parameters
---------------------------

This configuration submodule contains parameters related to the ortholog map pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{wildcards.sample}.genes.tab"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "ortholog_map"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.orthologs.tab"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {
    "memory": "40G"  # Custom memory request for ortholog mapping
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

resource_params = {}
"""
Dictionary defining the pipeline step's parameters that control resource usage but do not affect the output
"""

operating_params = {
    "type": "default",  # ID for operation to perform
    "method": "fractional"  # Ortholog counting method, defines how gene counts are distributed among orthologs they are annotated with, options include: fractional, whole
}
"""
Dictionary defining the pipeline step's parameters using when running the associated software
"""

benchmark_file = "{sample}.log"
"""
The benchmark filename pattern
"""

# Defining options for different operations to run during this step


def default(inputs, outputs, wildcards):
    """
    Default gene map operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    # Set locations of reference files if necessary
    gene_to_ortholog = fo.gene_to_ortholog_directory + op.gene_to_ortholog

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        subprocess.run([op.python, "src/ortholog_map.py", inputs.input, operating_params["method"], gene_to_ortholog, "--output", outputs[0]])

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
