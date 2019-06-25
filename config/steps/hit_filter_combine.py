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
    "options": "-cwd -R y -l disk_free=250G"  # Because combining read mapping output is not RAM intensive, but still generates large output before compression, reserve disk space in the job request so that we limit the number of hit filtering jobs based on available disk space
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

resource_params = {}
"""
Dictionary defining the pipeline step's parameters that control resource usage but do not affect the output
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

# Defining options for different operations to run during this step


def default(inputs, outputs, wildcards):
    """
    Default FASTQ summary operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    lf.combine_list_rows(inputs, outputs[0])


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
