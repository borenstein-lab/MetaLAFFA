from file_handling_lib import *

NUM_LINES_PER_READ = 4
READ_NAME_LINE = 0
QUALITY_LINE = 3


def get_base_quality_vector(quality_string):
    """
    Translates a string representing the quality of each base of a read into a numerical vector of quality scores

    :param quality_string: The string representing the quality of each base
    :return: A vector of numerical quality scores
    """

    # Subtract 33 to get the correct score
    return [x - 33 for x in map(ord, quality_string)]


def get_fastq_stats(fastq_filename):
    """
    Reads a FASTQ file and reports the number of reads, average read length, and average base quality

    :param fastq_filename: FASTQ file for which to generate summary statistics
    :return: A list of FASTQ summary statistics
    """

    # Initialize tracking variables to help parse the file
    line_count = 0
    read_count = 0
    base_count = 0
    total_quality_score = 0

    # Open a connection to the file
    f = custom_read(fastq_filename)

    # Iterate through the lines of the file
    for line in f:

        # If we are on a read name line, we found another read
        if line_count == READ_NAME_LINE:
            read_count += 1

        # Otherwise, if we are on the quality string line, we can get the length of the read and sum of base quality scores
        elif line_count == QUALITY_LINE:
            base_count += len(line.strip())
            total_quality_score += sum(get_base_quality_vector(line.strip()))

        # Increment our current line and reset if we've read all of the lines for a single read
        line_count += 1
        line_count = line_count % NUM_LINES_PER_READ
    f.close()

    return read_count, float(base_count) / float(read_count), float(total_quality_score) / float(base_count)
