"""
Duplicate filter step parameters
---------------------------

This configuration submodule contains parameters related to the duplicate filter pipeline step.
"""

import config.operation as op
import config.library_functions as lf
import config.file_organization as fo
import subprocess

input_dic = {
    "forward": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["forward"]},
    "reverse": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["reverse"]},
    "singleton": {"prior_step": 0, "file": "{wildcards.sample}.%s.fastq" % op.fastq_types["singleton"]}
}
"""
Dictionary defining the pipeline step's input structure.
"""

step_prefix = "duplicate_filter"
"""
The prefix to use in output subdirectory naming and provenance naming
"""

output_list = [
    "{sample}.%s.fastq" % op.fastq_types["forward"],
    "{sample}.%s.fastq" % op.fastq_types["reverse"],
    "{sample}.%s.fastq" % op.fastq_types["singleton"],
    "{sample}.R.marked",
    "{sample}.S.marked",
    "{sample}.R.metrics",
    "{sample}.S.metrics"
]
"""
List defining the pipeline step's outputs.
"""

cluster_params = {
    "memory": "40G"  # Custom memory request for duplicate filtering
}
"""
Dictionary defining the pipeline step's cluster parameters
"""

required_programs = {
    "picard": "picard",  # Path to PICARD tools executable
    "samtools": "samtools"  # Path to samtools executable
}
"""
Dictionary defining the paths to programs used by this pipeline step
"""

non_essential_params = {
    "mark_duplicates": "MarkDuplicates",  # PICARD tool for marking duplicates
    "fastq_to_sam": "FastqToSam",  # PICARD tool for converting FASTQ files to SAM format
    "quality_format": "Standard",  # SAM file quality format
    "sort_order": "coordinate"  # SAM file sort order
}
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

    # If either of the paired read files are non-empty, filter them for duplicate reads
    if not lf.is_empty(inputs.forward) or not lf.is_empty(inputs.reverse):

        # Convert the paired read files to SAM format
        paired_sam = inputs.forward + ".paired.sam"
        subprocess.run([required_programs["picard"], non_essential_params["fastq_to_sam"], "F1=%s" % inputs.forward, "F2=%s" % inputs.reverse, "O=%s" % paired_sam, "V=%s" % non_essential_params["quality_format"], "SO=%s" % non_essential_params["sort_order"], "SM=%s" % wildcards.sample])

        # Mark duplicates
        marked_output = inputs.forward + ".marked.sam"
        subprocess.run([required_programs["picard"], non_essential_params["mark_duplicates"], "I=%s" % paired_sam, "O=%s" % marked_output, "M=%s" % outputs[5]])
        subprocess.run(["rm", paired_sam])

        # Convert the marked output to a parse-able format
        formatted_marked_output = inputs.forward + ".samview"
        with open(formatted_marked_output, "w") as formatted_marked_output_file:
            subprocess.run([required_programs["samtools"], "view", marked_output], stdout=formatted_marked_output_file)
        subprocess.run(["rm", marked_output])

        # Extract duplicate reads from the formatted marked output
        with open(outputs[3], "w") as marked_read_file:
            subprocess.run([fo.source_directory + "extract_duplicates.py", formatted_marked_output], stdout=marked_read_file)

        # Remove the marked reads from the original FASTQs
        with open(outputs[0], "w") as forward_output:
            subprocess.run([fo.source_directory + "remove_marked_reads.py", outputs[3], inputs.forward], stdout=forward_output)
        with open(outputs[1], "w") as reverse_output:
            subprocess.run([fo.source_directory + "remove_marked_reads.py", outputs[3], inputs.reverse], stdout=reverse_output)

    # Otherwise, if both paired read files were dummy files, create dummy outputs
    else:
        subprocess.run(["touch", outputs[0], outputs[1], outputs[3], outputs[5]])

    # If the singleton read file is non-empty, filter it for duplicate reads
    if not lf.is_empty(inputs.singleton):

        # Convert the singleton read file to SAM format
        singleton_sam = inputs.forward + ".singleton.sam"
        subprocess.run([required_programs["picard"], non_essential_params["fastq_to_sam"], "F1=%s" % inputs.forward, "O=%s" % singleton_sam, "V=%s" % non_essential_params["quality_format"], "SO=%s" % non_essential_params["sort_order"], "SM=%s" % wildcards.sample])

        # Mark duplicates
        marked_output = inputs.singleton + ".marked.sam"
        subprocess.run([required_programs["picard"], non_essential_params["mark_duplicates"], "I=%s" % singleton_sam, "O=%s" % marked_output, "M=%s" % outputs[6]])
        subprocess.run(["rm", singleton_sam])

        # Convert the marked output to a parse-able format
        formatted_marked_output = inputs.singleton + ".samview"
        with open(formatted_marked_output, "w") as formatted_marked_output_file:
            subprocess.run([required_programs["samtools"], "view", marked_output], stdout=formatted_marked_output_file)
        subprocess.run(["rm", marked_output])

        # Extract duplicate reads from the formatted marked output
        with open(outputs[4], "w") as marked_read_file:
            subprocess.run([fo.source_directory + "extract_duplicates.py", formatted_marked_output], stdout=marked_read_file)

        # Remove the marked reads from the original FASTQs
        with open(outputs[2], "w") as singleton_output:
            subprocess.run([fo.source_directory + "remove_marked_reads.py", outputs[4], inputs.singleton], stdout=singleton_output)

    # Otherwise, if the singleton read file was a dummy file, create dummy outputs
    else:
        subprocess.run(["touch", outputs[2], outputs[4], outputs[6]])


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
