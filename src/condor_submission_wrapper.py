import argparse
import os
import re
import subprocess
import config.file_organization as fo
from snakemake.utils import read_job_properties

# Parse command line arguments
parser = argparse.ArgumentParser(description="Submits a Snakemake job to Condor")
parser.add_argument("job_script", help="Job script created by Snakemake")
args = parser.parse_args()

job_properties = read_job_properties(args.job_script)

# Create a directory to put Condor submission files if one does not already exist
if not os.path.isdir(fo.submission_file_directory):
    os.makedirs(fo.submission_file_directory)

cluster_params = job_properties["params"]
memory = cluster_params["memory"]
cores = cluster_params["cores"]

# Fix memory requirement to be in KB for Image_Size request
memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
if memory_units == "M":
    memory_int *= 10 ** 3
elif memory_units == "G":
    memory_int *= 10 ** 6

# Create a Condor submission file and submit the job to the cluster
submission_file_path = fo.submission_file_directory + os.path.basename(args.job_script)
with open(submission_file_path, "w") as submission_file:
    submission_file.write("Executable = %s\n" % args.job_script)
    # submission_file.write("Request_CPUs = %d\n" % cores)
    submission_file.write("Image_Size = %d\n" % memory_int)
    submission_file.write("Queue\n")

subprocess.run(["condor_submit", submission_file_path])
