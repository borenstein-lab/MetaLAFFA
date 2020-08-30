import os
import config.file_organization as fo

"""
Cluster parameters
---------------------

This configuration submodule contains definitions of how to generate cluster job requests and default parameters for resource requests.
"""

default_cluster_params = {
    "memory": "10G",
    "time": "24:0:0",
    "cores": 1,
    "wd": os.getcwd(),
    "reserve": True
}
"""
Dictionary defining the default pipeline step cluster parameters
"""

submission_wrapper = "src/sge_submission_wrapper.py"
"""
The path to the submission wrapper script for interpreting Snakemake-determined job parameters when using a cluster
"""

jobscript = "src/configured_jobscript.sh"
"""
The path to the jobscript that handles setup of parameters for cluster job session
"""
