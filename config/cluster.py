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


def sge_cluster_param_string_generator(memory, time, cpus, options):
    """
    Defines how to generate a cluster job request on an SGE cluster.

    :param memory: Amount of memory to request per core (format: <number><unit>, e.g. 10G for ten gigabytes of RAM)
    :param time: Amount of time to request for the job to run (format: <hours>:<minutes>:<seconds>, e.g. 24:0:0 for one day))
    :param cpus: Number of CPUs to request for the job to run (format: <number>, e.g. 2 for 2 CPUs)
    :param options: Other available options
    :return: A string that indicates the cluster resource requests
    """

    # If multiple CPUs are requested, then we need to request memory per CPU, not total memory
    if cpus > 1:
        memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
        memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
        memory_per_cpu = math.floor(float(memory_int) / float(cpus))
        memory = str(memory_per_cpu) + memory_units

    # Create the base resource request and add on the request for multiple cpus if necessary
    request = "%s %s=%s %s=%s" % (options, sge_memory_request_option, memory, sge_time_request_option, time)
    if cpus > 1:
        request += " -pe serial %d" % cpus

    return request


default_memory = "10G"
"""
Default amount of RAM required for individual pipeline cluster jobs.
"""

default_time = "24:0:0"
"""
Default time limit for individual pipeline cluster jobs.
"""

default_cpus = 1
"""
Default number of CPUs to request for individual pipeline cluster jobs.
"""

default_options = "-cwd -R y"
"""
Default options for individual pipeline cluster jobs.
"""
