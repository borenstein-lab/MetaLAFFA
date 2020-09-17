"""
Ortholog aggregation summary step parameters
---------------------------

This configuration submodule contains parameters related to the ortholog aggregation summary pipeline step.
"""

import config.file_organization as fo
import config.library_functions as lf
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{wildcards.mapping}.aggregated_orthologs.tab"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "ortholog_aggregation_summary"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{mapping}.summary.txt",
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
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
    "type": "default"  # ID for operation to perform
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
    Default gene map summary operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    if not lf.is_empty(inputs.input):
        subprocess.run([fo.source_directory + "ortholog_aggregation_summary.py",
                        inputs.input,
                        "--grouping_name", wildcards.mapping,
                        "--output", outputs[0]])
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
