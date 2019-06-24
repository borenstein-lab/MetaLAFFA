"""
Quality filter step parameters
---------------------------

This configuration submodule contains parameters related to the quality filter pipeline step.
"""

import config.operation as op
import config.library_functions as lf
import subprocess

input_dic = {
    "forward": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["forward"]},
    "reverse": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["reverse"]},
    "singleton": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["singleton"]}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "quality_filter"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.%s.fastq" % op.fastq_types["forward"],
    "{sample}.%s.fastq" % op.fastq_types["reverse"],
    "{sample}.%s.fastq" % op.fastq_types["singleton"],
    "{sample}.S_new.fastq",
    "{sample}.S_old.fastq"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
"""
Dictionary defining the pipeline step's cluster parameters
"""

resource_params = {
    "trimmer": "src/Trimmomatic-0.39/trimmomatic-0.39.jar",  # Location of the quality trimming and filtering program
}

operating_params = {
    "type": "default",  # ID for operation to perform
    "max_info": "90:0.7",  # Maximum info specification for trimmomatic
    "min_len": "60"  # Minimum read length to keep
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

    trimming_parameters = " ".join(["MAXINFO:" + operating_params["max_info"], "MINLEN:" + operating_params["min_len"]])

    # If either of the paired read files are non-empty, filter them for duplicate reads
    if not lf.is_empty(inputs.forward) or not lf.is_empty(inputs.reverse):

        # Assign intermediate files for the separate new singletons
        new_forward_singletons = outputs[3] + ".forward_singletons.fastq"
        new_reverse_singletons = outputs[3] + ".reverse_singletons.fastq"

        # Perform paired-end quality filtering and trimming
        subprocess.run([op.java, "-jar", resource_params["trimmer"], "PE", inputs.forward, inputs.reverse, outputs[0], new_forward_singletons, outputs[1], new_reverse_singletons, trimming_parameters])

        # Merge new singletons into single new singleton file
        with open(outputs[3], "w") as output_file:
            subprocess.run(["cat", new_forward_singletons, new_reverse_singletons], stdout=output_file)
        subprocess.run(["rm", new_forward_singletons, new_reverse_singletons])

        # Add new singletons to combined singleton output file
        with open(outputs[2], "w") as output_file:
            subprocess.run(["cat", outputs[3]], stdout=output_file)

    # Otherwise, if both paired read files were dummy files, create dummy outputs
    else:
        subprocess.run(["touch", outputs[0], outputs[1], outputs[3]])

    # If the singleton read file is non-empty, filter it for host reads
    if not lf.is_empty(inputs.singleton):

        # Perform single-end quality filtering and trimming
        subprocess.run([op.java, "-jar", resource_params["trimmer"], "SE", inputs.singleton, outputs[4],
                        trimming_parameters])

        # If we quality filtered the paired-end reads, add singletons to combined singleton output file
        if not lf.is_empty(inputs.forward) or not lf.is_empty(inputs.reverse):
            with open(outputs[2], "a") as output_file:
                subprocess.run(["cat", outputs[4]], stdout=output_file)

        # Otherwise, the combined singleton output file is the same as the results of the singleton quality filtering
        else:
            with open(outputs[2], "w") as output_file:
                subprocess.run(["cat", outputs[4]], stdout=output_file)

    # Otherwise, if the singleton read file was a dummy file, create dummy outputs
    else:
        subprocess.run(["touch", outputs[2], outputs[4]])


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
