# Author: Alex Eng
# Date: 11/21/2017
#
# A custom library for file handling functions

import gzip, re
from future import *

def custom_read(filename):
    '''Opens a file object appropriately depending on whether the file is compressed or not'''

    # If the file has been gzipped, open using the gzip library
    if re.search("\.gz$", filename):
        return gzip.open(filename, 'rt')

    # Otherwise, open normally
    else:
        return open(filename, 'r')
