"""
Ortholog aggregation step parameters
---------------------------

This configuration submodule contains parameters related to the ortholog aggregation pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "orthologs.tab"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "ortholog_aggregation"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{mapping}.aggregated_orthologs.tab"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
"""
Dictionary defining the pipeline step's cluster parameters
"""

resource_params = {}
"""
Dictionary defining the pipeline step's parameters that control resource usage but do not affect the output
"""

operating_params = {
    "type": "default",  # ID for operation to perform
    "method": "standard",  # Method to use for ortholog aggregation, options include: standard
    "standard_method": "fractional"  # Standard aggregation method to use, options include: fractional, whole
}
"""
Dictionary defining the pipeline step's parameters using when running the associated software
"""

benchmark_file = "{mapping}.log"
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

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        mapping = fo.ortholog_to_grouping_directory + wildcards.mapping


        if operating_params["method"] == "standard":
            subprocess.run([op.python, "src/ortholog_aggregation.py", inputs.input, operating_params["standard_method"], mapping, "--output", outputs[0]])

    # Otherwise, if the input file is a dummy file, create dummy outputs
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
