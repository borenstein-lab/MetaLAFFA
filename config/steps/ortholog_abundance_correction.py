"""
Ortholog abundance correction step parameters
---------------------------

This configuration submodule contains parameters related to the ortholog abundance correction pipeline step.
"""

from config import env
import config.library_functions as lf
import config.operation as op
import config.file_organization as fo
import subprocess

input_dic = {
    "input": {"prior_step": 0, "file": "orthologs.tab"}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "ortholog_abundance_correction"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "orthologs.tab"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
"""
Dictionary defining the pipeline step's cluster parameters
"""

required_programs = {
    "musicc": fo.python_source_directory + "run_musicc.py"  # Location of the MUSiCC program
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
    "method": "musicc",  # Method to use for ortholog abundance correction, options include: musicc, relative
    "musicc_method": ["-n", "-c", "learn_model"]  # MUSiCC correction method, options include: -n, -c use_generic, -c learn_model
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
    Default gene map operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    # If the input file is non-empty, map the reads
    if not lf.is_empty(inputs.input):

        if operating_params["method"] == "musicc":
            command = [required_programs["musicc"], inputs.input, "-o", outputs[0]] + operating_params["musicc_method"]
            subprocess.run(command, env=env)

        elif operating_params["method"] == "relative":
            subprocess.run([op.python, "src/ortholog_abundance_correction.py", inputs.input, "--output", outputs[0]], env=env)

        # Otherwise, if the method is unrecognized, just copy the input file to the output
        else:
            subprocess.run(["cp", inputs.input, outputs[0]], env=env)

    # Otherwise, if the input file is a dummy file, create a dummy output
    else:
        subprocess.run(["touch", outputs[0]], env=env)


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
