"""
Duplicate filter summary combine step parameters
---------------------------

This configuration submodule contains parameters related to the duplicate filter summary combine pipeline step.
"""

import config.file_organization as fo
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{sample}.{type}.summary.txt"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "duplicate_filter_summary_combine"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "summary.txt"
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

benchmark_file = "log"
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

    subprocess.run([fo.source_directory + "summary_combine.py"] + list(inputs.plainstrings()) + ["--output", outputs[0]])


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
