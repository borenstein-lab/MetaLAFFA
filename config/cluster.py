"""
Cluster parameters
---------------------

This configuration submodule contains parameters related to cluster interaction and resource requests.
"""

alignment_memory = 220
"""
Amount of RAM (in GB) required for each individual alignment cluster job. 
"""

alignment_cpus = 22
"""
Number of CPUs required for each individual alignment cluster job. 
"""

alignment_merging_disk_space = 250
"""
Amount of local disk space (in GB) required in the temporary file directory for merging alignment outputs. 
"""

alignment_filtering_disk_space = 200
"""
Amount of local disk space (in GB) required in the temporary file directory for filtering alignment output.
"""
