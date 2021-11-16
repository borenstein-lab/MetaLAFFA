"""
Hit filter combine step parameters
---------------------------

This configuration submodule contains parameters related to the hit filter combine pipeline step.
"""

import config.operation as op
import config.library_functions as lf

input_dic = {
    "forward": {"prior_step": 0, "file": "{wildcards.sample}.%s.mapping.txt" % op.fastq_types["forward"]},
    "reverse": {"prior_step": 0, "file": "{wildcards.sample}.%s.mapping.txt" % op.fastq_types["reverse"]},
    "singleton": {"prior_step": 0, "file": "{wildcards.sample}.%s.mapping.txt" % op.fastq_types["singleton"]}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "hit_filter_combine"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.mapping.txt"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {
    "disk": "250G"
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
    "type": "default"  # ID for operation to perform
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
    Default FASTQ summary operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :param log: The log file
    :return: None.
    """

    lf.combine_list_rows(list(inputs.plainstrings()), outputs[0], log[0])


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
