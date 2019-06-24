"""
Host filter step parameters
---------------------------

This configuration submodule contains parameters related to the host filter pipeline step.
"""

import config.operation as op
import config.file_organization as fo
import config.library_functions as lf
import os
import subprocess

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
    "{sample}.%s.fastq" % op.fastq_types["singleton"],
    "{sample}.R.marked.txt",
    "{sample}.S.marked.txt"
]
"""
List defining the pipeline step's output structure.
"""

cluster_params = {}
"""
Dictionary defining the pipeline step's cluster parameters
"""

resource_params = {
    "bmfilter": "src/bmtagger/bmtools/bmtagger/bmfilter",  # Path to bmfilter program
    "bmtagger": "src/bmtagger/bmtools/bmtagger/bmtagger.sh",  # Path to bmtagger program
    "extract_fa": "src/bmtagger/bmtools/bmtagger/extract_fullseq",  # Path to extract_fa program
    "srprism": "src/bmtagger/src/srprism/gnuac/app/srprism",  # Path to srprism program
    "blastn_dir": ":src/ncbi-blast-2.9.0+-src/c++/ReleaseMT/bin"  # Path to directory containing the blastn program
}
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

    # Set environment variables for bmtagger
    running_env = os.environ.copy()
    running_env["BMFILTER"] = resource_params["bmfilter"]
    running_env["EXTRACT_FA"] = resource_params["extract_fa"]
    running_env["SRPRISM"] = resource_params["srprism"]
    running_env["PATH"] += resource_params["blastn_dir"]

    # Set locations of reference files
    bitmask = fo.bitmask_directory + op.host_database + ".bitmask"
    host_database = fo.database_directory + op.host_database
    srprism_index = fo.srprism_index_directory + op.host_database + ".srprism"

    # If either of the paired read files are non-empty, filter them for host reads
    if not lf.is_empty(inputs.forward) or not lf.is_empty(inputs.reverse):

        # Match the read IDs between the paired read files for bmtagger
        matched_input = inputs.reverse + ".matched"
        with open(matched_input, "w") as matched_input_file:
            subprocess.run([op.python, "src/match_read_names.py", inputs.forward, inputs.reverse], stdout=matched_input_file)

        # bmtagger requires unzipped input files, so unzip the forward read file without deleting the original
        tmp_unzipped_exists = False
        forward_input = inputs.forward
        if lf.is_zipped(inputs.forward):
            with open(inputs.forward + "unzipped", "w") as tmp_unzipped_file:
                subprocess.run(["zcat", inputs.forward], stdout=tmp_unzipped_file)
            tmp_unzipped_exists = True
            forward_input = inputs.forward + "unzipped"

        # Run bmtagger on the paired (and read ID-matched) read files to mark reads
        subprocess.run([resource_params["bmtagger"], "-q1", "-1", forward_input, "-2", inputs.reverse + ".matched", "-o", outputs[3], "-b", bitmask, "-d", host_database, "-x", srprism_index], env=running_env)

        # Remove the temporary matched read ID file
        subprocess.run(["rm", matched_input])

        # Remove the temporary unzipped forward read file if it exists
        if tmp_unzipped_exists:
            subprocess.run(["rm", forward_input])

        # Remove the marked reads from the original fastqs
        with open(outputs[0], "w") as forward_output:
            subprocess.run([op.python, "src/remove_marked_reads.py", outputs[3], inputs.forward], stdout=forward_output)
        with open(outputs[1], "w") as reverse_output:
            subprocess.run([op.python, "src/remove_marked_reads.py", outputs[3], inputs.reverse],
                           stdout=reverse_output)

    # Otherwise, if both paired read files were dummy files, create dummy outputs
    else:
        subprocess.run(["touch", outputs[0], outputs[1], outputs[3]])

    # If the singleton read file is non-empty, filter it for host reads
    if not lf.is_empty(inputs.singleton):

        # bmtagger requires unzipped input files, so unzip the forward read file without deleting the original
        tmp_unzipped_exists = False
        singleton_input = inputs.singleton
        if lf.is_zipped(inputs.singleton):
            with open(inputs.singleton + "unzipped", "w") as tmp_unzipped_file:
                subprocess.run(["zcat", inputs.singleton], stdout=tmp_unzipped_file)
            tmp_unzipped_exists = True
            singleton_input = inputs.singleton + "unzipped"

        # Run bmtagger on the singleton read file to mark reads
        subprocess.run([resource_params["bmtagger"], "-q1", "-1", singleton_input, "-o", outputs[4], "-b", bitmask, "-d", host_database, "-x", srprism_index], env=running_env)

        # Remove the temporary unzipped forward read file if it exists
        if tmp_unzipped_exists:
            subprocess.run(["rm", singleton_input])

        # Remove the marked reads from the original fastqs
        with open(outputs[2], "w") as singleton_output:
            subprocess.run([op.python, "src/remove_marked_reads.py", outputs[4], inputs.singleton], stdout=singleton_output)

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
