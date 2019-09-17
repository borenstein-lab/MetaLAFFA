"""
Host filter step parameters
---------------------------

This configuration submodule contains parameters related to the host filter pipeline step.
"""

from config import env
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

cluster_params = {
    "memory": "20G"  # Custom memory requirement for host filtering
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

required_programs = {
    "bmtool": fo.source_directory + "bmtools/bmtagger/bmtool",  # Path to bmtools program
    "bmfilter": fo.source_directory + "bmtools/bmtagger/bmfilter",  # Path to bmfilter program
    "bmtagger": fo.source_directory + "bmtools/bmtagger/bmtagger.sh",  # Path to bmtagger program
    "extract_fa": fo.source_directory + "bmtools/bmtagger/extract_fullseq",  # Path to extract_fa program
    "srprism": fo.source_directory + "srprism/gnuac/app/srprism",  # Path to srprism program
    "blastn": fo.source_directory + "ncbi-blast-2.2.31+-src/c++/ReleaseMT/bin/blastn"  # Path to blastn program
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

    # Set environment variables for bmtagger
    running_env = env
    running_env["BMFILTER"] = required_programs["bmfilter"]
    running_env["EXTRACT_FA"] = required_programs["extract_fa"]
    running_env["SRPRISM"] = required_programs["srprism"]
    running_env["BLASTN"] = required_programs["blastn"]

    # Set locations of reference files
    bitmask = op.host_bitmask_file
    host_database = op.host_database_file
    srprism_index = op.host_index_file

    # If either of the paired read files are non-empty, filter them for host reads
    if not lf.is_empty(inputs.forward) or not lf.is_empty(inputs.reverse):

        # Match the read IDs between the paired read files for bmtagger
        matched_input = inputs.reverse + ".matched"
        with open(matched_input, "w") as matched_input_file:
            subprocess.run([op.python, fo.source_directory + "match_read_names.py", inputs.forward, inputs.reverse], stdout=matched_input_file, env=env)

        # bmtagger requires unzipped input files, so unzip the forward read file without deleting the original
        tmp_unzipped_exists = False
        forward_input = inputs.forward
        if lf.is_zipped(inputs.forward):
            with open(inputs.forward + "unzipped", "w") as tmp_unzipped_file:
                subprocess.run(["zcat", inputs.forward], stdout=tmp_unzipped_file, env=env)
            tmp_unzipped_exists = True
            forward_input = inputs.forward + "unzipped"

        # Run bmtagger on the paired (and read ID-matched) read files to mark reads
        subprocess.run([required_programs["bmtagger"], "-q1", "-1", forward_input, "-2", inputs.reverse + ".matched", "-o", outputs[3], "-b", bitmask, "-d", host_database, "-x", srprism_index], env=running_env)

        # Remove the temporary matched read ID file
        subprocess.run(["rm", matched_input], env=env)

        # Remove the temporary unzipped forward read file if it exists
        if tmp_unzipped_exists:
            subprocess.run(["rm", forward_input], env=env)

        # Remove the marked reads from the original fastqs
        with open(outputs[0], "w") as forward_output:
            subprocess.run([op.python, fo.source_directory + "remove_marked_reads.py", outputs[3], inputs.forward], stdout=forward_output, env=env)
        with open(outputs[1], "w") as reverse_output:
            subprocess.run([op.python, fo.source_directory + "remove_marked_reads.py", outputs[3], inputs.reverse], stdout=reverse_output, env=env)

    # Otherwise, if both paired read files were dummy files, create dummy outputs
    else:
        subprocess.run(["touch", outputs[0], outputs[1], outputs[3]], env=env)

    # If the singleton read file is non-empty, filter it for host reads
    if not lf.is_empty(inputs.singleton):

        # bmtagger requires unzipped input files, so unzip the forward read file without deleting the original
        tmp_unzipped_exists = False
        singleton_input = inputs.singleton
        if lf.is_zipped(inputs.singleton):
            with open(inputs.singleton + "unzipped", "w") as tmp_unzipped_file:
                subprocess.run(["zcat", inputs.singleton], stdout=tmp_unzipped_file, env=env)
            tmp_unzipped_exists = True
            singleton_input = inputs.singleton + "unzipped"

        # Run bmtagger on the singleton read file to mark reads
        subprocess.run([required_programs["bmtagger"], "-q1", "-1", singleton_input, "-o", outputs[4], "-b", bitmask, "-d", host_database, "-x", srprism_index], env=running_env)

        # Remove the temporary unzipped forward read file if it exists
        if tmp_unzipped_exists:
            subprocess.run(["rm", singleton_input], env=env)

        # Remove the marked reads from the original fastqs
        with open(outputs[2], "w") as singleton_output:
            subprocess.run([op.python, fo.source_directory + "remove_marked_reads.py", outputs[4], inputs.singleton], stdout=singleton_output, env=env)

    # Otherwise, if the singleton read file was a dummy file, create dummy outputs
    else:
        subprocess.run(["touch", outputs[2], outputs[4]], env=env)


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
