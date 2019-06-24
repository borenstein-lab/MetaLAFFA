"""
Library functions
---------------------

This configuration submodule consolidates functions used by all parts of the pipeline to standardize how operations are handled.
"""

import config.operation as op
import config.file_organization as fo
import re
import os
import subprocess
import sys
import copy


def unzipped_name(filename):
    """
    Determines the unzipped version of a file name (remove the ".gz" suffix if present).

    :param filename: File name to "unzip"
    :return: Unzippped file name
    """

    return re.sub("\\.gz$", "", filename)


def get_working_output_name(filename):
    """
    Determines the version of the file name that will be operated on (e.g. if operations should be performed in a temporary directory, rather than in the output directory).

    :param filename: File name to process
    :return: Working name of the file
    """

    working_output = unzipped_name(filename)
    if op.work_in_tmp_dir:
        working_output = fo.tmp_dir + os.path.basename(working_output)
    return working_output


def process_output(filename):
    """
    Standard output file processing for pipeline steps (e.g. if specified, zip output files).

    :param filename: File to process
    :return: None
    """

    working_output = get_working_output_name(filename)
    if op.zipped_files:
        subprocess.run(["gzip", working_output])
        working_output = working_output + ".gz"
    if op.work_in_tmp_dir:
        subprocess.run(["mv", working_output, filename])


def is_zipped(filename):
    """
    Tests if a file name represents a zipped file.

    :param filename: File name to test
    :return: Whether the file is zipped or not
    """

    return re.search("\\.gz$", filename) is not None


def is_empty(filename):
    """
    Returns whether a file is empty

    :param filename: File to check
    :return: True or False for whether the file is empty
    """

    num_lines = None

    # If the file is zipped, then we need to unzip it before checking the contents
    if is_zipped(filename):

        # Unzip the file and then count the number of lines
        unzip_process = subprocess.Popen(["zcat", filename], stdout=subprocess.PIPE)
        num_lines = int(subprocess.run(["wc -l"], stdin=unzip_process.stdout, capture_output=True, text=True).stdout.strip())

    # Otherwise, we just count the number of lines
    else:
        num_lines = int(subprocess.run(["wc -l", filename], capture_output=True, text=True).stdout.strip())

    # Return whether the number of lines is 0 (an empty file)
    return num_lines == 0


def combine_list_rows(lists_to_combine, combined_list):
    """
    Standard combining of list rows.

    :param lists_to_combine: List of list files to combine
    :param combined_list: Name of output list file
    :return: None
    """

    # Initialize the combined list file
    with open(combined_list, "w") as combined_list_file:
        if is_zipped(lists_to_combine[0]):
            subprocess.run(["zcat", lists_to_combine[0]], stdout=combined_list_file)
        else:
            subprocess.run(["cat", lists_to_combine[0]], stdout=combined_list_file)

    # Now add rows from the rest of the lists
    with open(combined_list, "a") as combined_list_file:
        for list in lists_to_combine[1:]:
            if is_zipped(list):
                subprocess.run(["zcat", list], stdout=combined_list_file)
            else:
                subprocess.run(["cat", list], stdout=combined_list_file)


def replace_restricted_patterns(pattern, pattern_restrictions):
    """
    Replaces wildcards in a pattern with a regex restricting that wildcard to a set of possible patterns

    :param pattern: Pattern to modify
    :param pattern_restrictions: Dictionary of wildcard restrictions
    :return: Modified pattern with restricted wildcards
    """

    for wildcard in pattern_restrictions:
        pattern = re.sub("{" + wildcard + "}", "(" + "|".join(pattern_restrictions[wildcard]) + ")", pattern)
    return pattern


def get_initial_data_patterns(step_info, step_params):
    """
    Generates a list of patterns used in input files.

    :param step_info: Dictionary describing which step outputs are used as inputs for other steps
    :param step_params: Dictionary describing the parameters of each step
    :return:
    """

    initial_steps = []
    for step_name in step_info:
        if "INPUT" in step_info[step_name]:
            initial_steps.append(step_name)
    initial_patterns = []
    for step_name in initial_steps:
        for input_name in step_params[step_name]["input"](None):
            initial_patterns.append(step_params[step_name]["input"](None)[input_name])
    return initial_patterns


def generate_sample_list(step_info, step_params):
    """
    Generate the list of samples that need to be processed by the pipeline.

    :param step_info: Dictionary describing which step outputs are used as inputs for other steps
    :param step_params: Dictionary describing the parameters of each step
    :return: List of samples to process
    """

    # Initialize the container of samples to process as a set for initial creation without duplication, later convert to list for iteration
    samples = set()

    # Get the initial data file patterns from pipeline steps that use initial data as input
    initial_patterns = get_initial_data_patterns(step_info, step_params)

    # If a sample list is specified, process the samples in the list rather than determining the list of samples from the input data directory
    if op.sample_list is not None:

        # Read sample IDs from the file and add any with existing input files to the list of samples to process
        with open(op.sample_list) as sample_list_file:
            for line in sample_list_file:
                sample = line.strip()

                # Check for sample ID matches against input file patterns
                input_exists = False
                for input_file in os.listdir(fo.initial_data_directory):
                    for initial_pattern in initial_patterns:
                        fixed_pattern = re.sub("wildcards\\.", "", initial_pattern)
                        fixed_pattern = re.sub("{sample}", sample, fixed_pattern)
                        fixed_pattern = re.sub("\\.", "\\.", fixed_pattern)
                        fixed_pattern = replace_restricted_patterns(fixed_pattern, op.wildcard_restrictions)
                        if re.match(fixed_pattern, fo.initial_data_directory + input_file):
                            input_exists = True
                            break
                    if input_exists:
                        break

                # If input files exist for this sample, add it to the list
                if input_exists:
                    samples.add(sample)

                # Otherwise, inform the user that this sample has no associated input files
                else:
                    print(sample + " has no associated input files.")

        # If there are no samples in the sample list, report this to the user and exit
        if (len(samples)) == 0:
            sys.exit("There were no input files matching the sample IDs specified in the sample list.")

    # Otherwise, determine the list of samples from those files
    else:

        # Check for input files that match the expected naming format and extract sample IDs from those file names
        for input_file in os.listdir(fo.initial_data_directory):
            for initial_pattern in initial_patterns:
                fixed_pattern = re.sub("wildcards\\.", "", initial_pattern)
                fixed_pattern = re.sub("\\.", "\\.", fixed_pattern)
                fixed_pattern = replace_restricted_patterns(fixed_pattern, op.wildcard_restrictions)
                fixed_pattern = re.sub("{sample}", "([^.]*)", fixed_pattern)
                match = re.match(fixed_pattern, fo.initial_data_directory + input_file)
                if match is not None:
                    samples.add(match.group(1))

        # If there are no samples in the sample list, report this to the user and exit
        if (len(samples)) == 0:
            sys.exit("There were no input files matching the expected input file name format.")

    # Convert the samples container to a list for iteration
    samples = list(samples)

    return samples


def generate_possible_patterns_from_restricted_wildcards(initial_patterns, wildcard_restriction_dic):
    """
    Takes patterns containing wildcards replaces wildcards that have been restricted to a set of specific set of values to generate all possible sub-patterns that use those restricted wildcard values.

    :param initial_patterns: The template patterns used to generate the possible sub-patterns
    :param wildcard_restriction_dic: A dictionary that gives the restricted values for one or more wildcards
    :return: A list of patterns using all possible combinations of restricted wildcard values
    """

    possible_patterns = []
    restricted_wildcards = list(wildcard_restriction_dic.keys())

    # Generate sub-patterns for each provided initial pattern
    for initial_pattern in initial_patterns:

        # Remove the "wildcards" prefix, which is only used for Snakemake unpacking
        fixed_pattern = re.sub("wildcards\\.", "", initial_pattern)

        # We need to fill in all of the possible wild card combinations from the restricted wild card list
        wildcard_value_indices = [0] * len(restricted_wildcards)
        while wildcard_value_indices[0] > -1:

            # Create a dictionary to format the pattern with the wildcard values
            formatting_dic = {}
            for wildcard_index in range(len(restricted_wildcards)):
                wildcard = restricted_wildcards[wildcard_index]
                formatting_dic[wildcard] = wildcard_restriction_dic[wildcard][wildcard_value_indices[wildcard_index]]

            # Format the pattern and add it to the output list
            possible_patterns.append(fixed_pattern.format(**formatting_dic))

            # Increment the wildcard value indices
            increment_index = len(wildcard_value_indices) - 1
            incremented = False
            while not incremented and increment_index >= 0:

                # If the current wildcard value index is not at the last value for this wildcard, increment this wildcard value index and we're done
                if wildcard_value_indices[increment_index] < len(wildcard_restriction_dic[restricted_wildcards[increment_index]]) - 1:
                    wildcard_value_indices[increment_index] += 1
                    incremented = True

                # Otherwise, roll it over to the first value for this wildcard and move on to the previous wildcard value index
                else:
                    wildcard_value_indices[increment_index] = 0
                    increment_index -= 1

            # If we didn't increment, then we've gone through all of the possible restricted wildcard combinations and we're done
            if not incremented:
                wildcard_value_indices[0] = -1

    return possible_patterns


def create_dummy_inputs(sample_list, step_info, step_params):
    """
    Creates dummy input files if necessary

    :param sample_list: List of samples
    :param step_info: Dictionary describing which step outputs are used as inputs for other steps
    :param step_params: Dictionary describing the parameters of each step
    :return: None
    """

    # Get initial data file patterns
    initial_patterns = get_initial_data_patterns(step_info, step_params)

    # Generate the expected input file names
    wildcard_restrictions_with_samples = copy.deepcopy(op.wildcard_restrictions)
    wildcard_restrictions_with_samples["sample"] = sample_list

    expected_input_files = generate_possible_patterns_from_restricted_wildcards(initial_patterns, wildcard_restrictions_with_samples)

    # For each sample, if the initial pattern includes a "{sample}" wildcard, create dummy files for any missing initial patterns and write which dummy files were created in a log file
    with open(op.dummy_input_log, "a") as dummy_log_file:

        # Check whether each expected input file exists, and if not create a dummy file and note that in the log file
        for expected_input_file in expected_input_files:
            if not os.path.isfile(expected_input_file):
                subprocess.run(["touch", get_working_output_name(expected_input_file)])
                process_output(expected_input_file)
                dummy_log_file.write(expected_input_file + os.linesep)


def run_step(curr_step_params, inputs, outputs, wildcards):
    """
    Runs standard and rule-specific operations

    :param curr_step_params: Dictionary of parameters for the current pipeline step
    :param inputs: The Snakemake inputs object
    :param outputs: The dictionary of output files
    :param wildcards: The Snakemake wildcards object
    :return: None
    """

    working_outputs = []
    for output in outputs:
        working_outputs.append(get_working_output_name(output))
    curr_step_params["rule_function"](inputs, working_outputs, wildcards)
    for output in outputs:
        process_output(output)
