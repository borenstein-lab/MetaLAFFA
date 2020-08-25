"""
Host filter step parameters
---------------------------

This configuration submodule contains parameters related to the host filter pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import subprocess
import re

input_dic = {
    "forward": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["forward"]},
    "reverse": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["reverse"]},
    "singleton": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["singleton"]}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "host_filter"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.%s.fastq" % op.fastq_types["forward"],
    "{sample}.%s.fastq" % op.fastq_types["reverse"],
    "{sample}.%s.fastq" % op.fastq_types["singleton"]
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {
    "memory": "20G"  # Custom memory requirement for host filtering
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

required_programs = {
    "bowtie2": "bowtie2",  # Path to Bowtie 2 program
}
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

# Defining options for different operations to run during this step


def default(inputs, outputs, wildcards):
    """
    Default FASTQ summary operations.

    :param inputs: Object containing the input file names
    :param outputs: Dictionary containing the output file names
    :param wildcards: Wildcards determined from input file name patterns
    :return: None.
    """

    # Set locations of reference files
    host_index = op.host_index

    # If either of the paired read files are non-empty, filter them for host reads
    if not lf.is_empty(inputs.forward) or not lf.is_empty(inputs.reverse):

        # Run Bowtie 2 on the paired read files
        subprocess.run([required_programs["bowtie2"],
                        "-x", host_index,
                        "-1", inputs.forward,
                        "-2", inputs.reverse,
                        "--un-conc", outputs[0],
                        "--no-unal",
                        "-S", "/dev/null"])

        # Fix output file names
        subprocess.run(["mv", re.sub('(\\.[^.]*)$', r'.1\1', outputs[0]), outputs[0]])
        subprocess.run(["mv", re.sub('(\\.[^.]*)$', r'.2\1', outputs[0]), outputs[1]])

    # Otherwise, if both paired read files were dummy files, create dummy outputs
    else:
        subprocess.run(["touch", outputs[0], outputs[1]])

    # If the singleton read file is non-empty, filter it for host reads
    if not lf.is_empty(inputs.singleton):

        # Run Bowtie 2 on the singleton read file
        subprocess.run([required_programs["bowtie2"],
                        "-x", host_index,
                        "-U", inputs.singleton,
                        "--un", outputs[2],
                        "--no-unal",
                        "-S", "/dev/null"])

    # Otherwise, if the singleton read file was a dummy file, create dummy outputs
    else:
        subprocess.run(["touch", outputs[2]])


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
