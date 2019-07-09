"""
Cluster parameters
---------------------

This configuration submodule contains definitions of how to generate cluster job requests and default parameters for resource requests.
"""

import re
import math

sge_memory_request_option = "-l mfree"
"""
The option used in SGE cluster job submissions to request RAM.
"""

sge_time_request_option = "-l h_rt"
"""
The option used in SGE cluster job submissions to request time.
"""


def sge_cluster_param_string_generator(memory, time, cores, options):
    """
    Defines how to generate a cluster job request on an SGE cluster.

    :param memory: Amount of memory to request per core (format: <number><unit>, e.g. 10G for ten gigabytes of RAM)
    :param time: Amount of time to request for the job to run (format: <hours>:<minutes>:<seconds>, e.g. 24:0:0 for one day))
    :param cores: Number of cores to request for the job to run (format: <number>, e.g. 2 for 2 cores)
    :param options: Other available options
    :return: A string that indicates the cluster resource requests
    """

    # If multiple cores are requested, then we need to request memory per CPU, not total memory
    if cores > 1:
        memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
        memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
        memory_per_cpu = math.floor(float(memory_int) / float(cores))
        memory = str(memory_per_cpu) + memory_units

    # Create the base resource request and add on the request for multiple cores if necessary
    request = "%s %s=%s %s=%s" % (options, sge_memory_request_option, memory, sge_time_request_option, time)
    if cores > 1:
        request += " -pe serial %d" % cores

    return request


default_memory = "10G"
"""
Default amount of RAM required for individual pipeline cluster jobs.
"""

default_time = "24:0:0"
"""
Default time limit for individual pipeline cluster jobs.
"""

default_cores = 1
"""
Default number of cores to request for individual pipeline cluster jobs.
"""

default_options = "-cwd -R y"
"""
Default options for individual pipeline cluster jobs.
"""
