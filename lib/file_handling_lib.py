import gzip
import re


def custom_read(filename):
    """
    Opens a file object appropriately depending on whether the file is compressed or not

    :param filename: File to open.
    :return: File object.
    """

    # If the file has been gzipped, open using the gzip library
    if re.search("\\.gz$", filename):
        return gzip.open(filename, 'rt')

    # Otherwise, open normally
    else:
        return open(filename, 'r')
