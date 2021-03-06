#!/usr/bin/env python

import argparse
import re
import math
import subprocess
from snakemake.utils import read_job_properties

# Parse command line arguments
parser = argparse.ArgumentParser(description="Submits a Snakemake job to SGE")
parser.add_argument("job_script", help="Job script created by Snakemake")
args = parser.parse_args()

job_properties = read_job_properties(args.job_script)

# Create an SGE submission command and submit the job to the cluster
cluster_params = job_properties["params"]["cluster"]
memory = cluster_params["memory"]
time = cluster_params["time"]
cores = cluster_params["cores"]
wd = cluster_params["wd"]
reserve = cluster_params["reserve"]
disk = None
if "disk" in cluster_params:
    disk = cluster_params["disk"]

#############################################
# EDIT BELOW HERE FOR CLUSTER CUSTOMIZATION #
#############################################

submission_command = ["qsub"]

# If multiple cores are requested, then we need to request memory per CPU, not total memory
if cores > 1:
    memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
    memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
    memory_per_cpu = math.floor(float(memory_int) / float(cores))
    memory = str(memory_per_cpu) + memory_units
    submission_command += ["-pe", "serial", str(cores)]

submission_command += ["-l", "mfree=%s" % memory, "-l", "h_rt=%s" % time, "-wd", wd]

# Only add the disk_free request if working in local TMP storage, otherwise comment out this addition to the resource request
if disk is not None:
    submission_command += ["-l", "disk_free=%s" % disk]

# Should resources for this job be reserved (set aside as they become available from completed jobs) until enough resources are available for the job to run?
if reserve:
    submission_command += ["-R", "y"]

submission_command += [args.job_script]

subprocess.run(submission_command)
