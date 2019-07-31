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

submission_command = ["qsub", args.job_script]

# If multiple cores are requested, then we need to request memory per CPU, not total memory
if cores > 1:
    memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
    memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
    memory_per_cpu = math.floor(float(memory_int) / float(cores))
    memory = str(memory_per_cpu) + memory_units
    submission_command += ["-pe serial %d" % cores]

submission_command += ["-l mfree=%s" % memory, "-l h_rt=%s" % time, "-wd %s" % wd]

if reserve:
    submission_command += ["-R y"]

subprocess.run(submission_command)
