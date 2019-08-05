"""
Summary combine step parameters
---------------------------

This configuration submodule contains parameters related to the summary combine pipeline step.
"""

import config.file_organization as fo
import config.library_functions as lf
import re
import os
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "FINAL_OUTPUTS"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "process_final_output"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = ["FINAL_OUTPUTS"]
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

benchmark_file = "log"
"""
The benchmark filename pattern
"""

# Defining options for different operations to run during this step


def default(inputs, outputs, wildcards):
    """
    Default final output processing operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    for final_output in list(inputs.plainstrings()):
        subprocess.run(["cp", final_output, lf.process_final_output_name(final_output)])


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
