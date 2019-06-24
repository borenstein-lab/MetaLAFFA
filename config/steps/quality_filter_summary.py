"""
Quality filter summary step parameters
---------------------------

This configuration submodule contains parameters related to the quality filter summary pipeline step.
"""

import config.operation as op
import subprocess

input_dic = {
    "pre_forward": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["forward"]},
    "pre_reverse": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["reverse"]},
    "post_forward": {"prior_step": 1, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["forward"]},
    "post_reverse": {"prior_step": 1, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["reverse"]},
    "new_singleton": {"prior_step": 1, "file": "{wildcards.sample}.S_new.fastq"},
    "old_singleton": {"prior_step": 1, "file": "{wildcards.sample}.S_old.fastq"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "quality_filter_summary"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.summary.txt",
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
"""
Dictionary defining the pipeline step's cluster parameters
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

    subprocess.run([op.python, "src/quality_filter_summary.py", inputs.pre_forward, inputs.pre_reverse, inputs.post_forward, inputs.post_reverse, inputs.new_singleton, inputs.old_singleton, "--output", outputs[0], "--use-sample"])


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
