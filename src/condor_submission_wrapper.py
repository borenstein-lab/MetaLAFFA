import argparse
import os
import re
import subprocess
from snakemake.utils import read_job_properties
from config import env

# Parse command line arguments
parser = argparse.ArgumentParser(description="Submits a Snakemake job to Condor")
parser.add_argument("job_script", help="Job script created by Snakemake")
args = parser.parse_args()

job_properties = read_job_properties(args.job_script)

# Create a directory to put Condor submission files if one does not already exist
submission_dir = "submission_files/"
if not os.path.isdir(submission_dir):
    os.makedirs(submission_dir)

cluster_params = job_properties["params"]["cluster"]
memory = cluster_params["memory"]
cores = cluster_params["cores"]

#############################################
# EDIT BELOW HERE FOR CLUSTER CUSTOMIZATION #
#############################################

# If your MetaLAFFA installation is not going to be available on all cluster nodes (e.g. you own a subset of nodes that, by default, mount the drive where you installed MetaLAFFA), then you will need to restrict Condor to requesting slots on the specific node or nodes that will have access to MetaLAFFA. If a single node, set node to the name of the node (in string format). If multiple nodes, set node to the OR'd union (using '|' to join) of those node names (in string format). For example, if you want to request slots on either "node1" or "node2", you would set: node = "node1 | node2".
node = None

# Fix memory requirement to be in KB for Image_Size request
memory_int = int(re.match("([0-9]*)[^0-9]", memory).group(1))
memory_units = re.match("[0-9]*([^0-9]*)$", memory).group(1)
if memory_units == "M":
    memory_int *= 10 ** 3
elif memory_units == "G":
    memory_int *= 10 ** 6

# Create a Condor submission file and submit the job to the cluster)
submission_file_path = submission_dir + os.path.basename(args.job_script)
with open(submission_file_path, "w") as submission_file:
    submission_file.write("Executable = %s\n" % args.job_script)
    if node is not None:
        submission_file.write("Requirements = (Machine == \"%s\")\n" % node)
    # submission_file.write("Request_CPUs = %d\n" % cores)
    submission_file.write("Image_Size = %d\n" % memory_int)
    submission_file.write("Queue\n")

subprocess.run(["condor_submit", submission_file_path], env=env)
