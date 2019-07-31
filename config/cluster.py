"""
Cluster parameters
---------------------

This configuration submodule contains definitions of how to generate cluster job requests and default parameters for resource requests.
"""

default_cluster_params = {
    "memory": "10G",
    "time": "24:0:0",
    "cores": 1,
    "wd": ".",
    "reserve": True
}
"""
Dictionary defining the default pipeline step cluster parameters
"""