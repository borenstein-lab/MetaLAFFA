"""
Hit filter step parameters
---------------------------

This configuration submodule contains parameters related to the hit filter pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "{wildcards.sample}.{wildcards.type}.mapping.txt"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "hit_filter"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.{type}.mapping.txt"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {
    "disk": "200G"  # Because hit filtering is not RAM intensive, but still generates large output before compression, reserve disk space in the job request so that we limit the number of hit filtering jobs based on available disk space
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
    "type": "default",  # ID for operation to perform
    "method": "best_hit",  # Hit filtering method (if using default hit filtering program, options include: best_hit, best_ortholog, best_n_hits, best_n_orthologs
    # "best_n": 10  # If keeping the top N, rather than just best, hits, this is how many to keep (comment out if not using a top N method, otherwise this will show up in the provenance as one of the settings)

}
"""
Dictionary defining the pipeline step's parameters using when running the associated software
"""

benchmark_file = "{sample}.{type}.log"
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

    # Set locations of reference files if necessary
    n_method = operating_params["method"] in ["best_n_hits", "best_n_orthologs"]
    ortholog_method = operating_params["method"] in ["best_ortholog", "best_n_orthologs"]
    gene_to_ortholog = None
    if ortholog_method:
        gene_to_ortholog = op.gene_to_ortholog_file

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        # Create the shell command to run, adding optional parameters as necessary
        command = [op.python, "src/hit_filter.py", inputs.input, operating_params["method"], "--output", outputs[0]]
        if ortholog_method:
            command += ["--gene_to_ortholog_map", gene_to_ortholog]
        if n_method:
            command += ["-n", operating_params["best_n"]]

        subprocess.run(command)

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
