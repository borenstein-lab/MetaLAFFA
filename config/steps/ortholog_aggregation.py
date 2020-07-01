"""
Ortholog aggregation step parameters
---------------------------

This configuration submodule contains parameters related to the ortholog aggregation pipeline step.
"""

from config import env
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
    "{mapping}.aggregated_orthologs.tab",
    "{mapping}.mapped_ortholog_counts.tab"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
"""
Dictionary defining the pipeline step's cluster parameters
"""

required_programs = {
    "empanada": fo.python_source_directory + "run_empanada.py"  # Location of the EMPANADA program
}
"""
Dictionary defining the paths to programs used by this pipeline step
"""

non_essential_params = {}
"""
Dictionary defining the pipeline step's parameters that don't affect the output
"""

operating_params = {
    "type": "default",  # ID for operation to perform
    "method": "empanada",  # Method to use for ortholog aggregation, options include: empanada, fractional
    "empanada_method": ["-map", "by_support"]  # EMPANADA aggregation method, options include: -map by_support, "
                                          "-map naive, -map by_sum_abundance, -map by_avg_abundance
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
    Default ortholog aggregation operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        mapping = fo.ortholog_to_grouping_directory + wildcards.mapping + op.ortholog_to_grouping_suffix

        if operating_params["method"] == "empanada":
            command = [required_programs["empanada"], "-ko", inputs.input, "-ko2path", mapping, "-o", outputs[0], "-oc", outputs[1]] + operating_params["empanada_method"]
            subprocess.run(command, env=env)

        # Otherwise, if the method is unrecognized, just copy the input file to the output
        else:
            subprocess.run(["cp", inputs.input, outputs[0]], env=env)
            subprocess.run(["touch", outputs[1]], env=env)

    # Otherwise, if the input file is a dummy file, create dummy outputs
    else:
        subprocess.run(["touch", outputs[0], outputs[1]], env=env)


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
