"""
Gene map summary step parameters
---------------------------

This configuration submodule contains parameters related to the gene map summary pipeline step.
"""

import config.operation as op
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{wildcards.sample}.genes.tab"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "gene_map_summary"
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
    Default gene map summary operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    subprocess.run([op.python, "src/gene_map_summary.py", inputs.input, "--output", outputs[0]])


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
