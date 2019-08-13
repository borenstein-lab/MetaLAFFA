import argparse
from config import env
import subprocess
import config.operation as op

parser = argparse.ArgumentParser(description="Wrapper script to run the pipeline, ensuring that the shell environment for Snakemake is correct.")
parser.add_argument("--submission_wrapper", "-s", default="python3 src/sge_submission_wrapper.py", help="The wrapper script to use to parse job settings and resource requests and submit a job to the cluster (default %(default)s)")
parser.add_argument("--number_of_jobs", "-n", type=int, help="Number of jobs have submitted to the cluster queue at any one time, or if running locally, the number of cores to use for running jobs in parallel.")
parser.add_argument("--wait", "-w", type=int, help="Number of seconds to wait after a job finishes before checking that the output exists, this can avoid Snakemake incorrectly marking a step as failed when a file might not be present due to latency over a shared file system.")
parser.add_argument("--local", "-l", action="store_true", help="If used, run metaLAFFA locally rather than on a cluster.")

args = parser.parse_args()

# Create the snakemake command and add options to it as specified by the arguments
command = [op.snakemake]
if not args.local:
    command += ["-c", args.submission_wrapper]
if args.number_of_jobs is not None:
    command += ["-j", str(args.number_of_jobs)]
if args.wait is not None:
    command += ["--latency-wait", str(args.wait)]

# Run the snakemake command
subprocess.run(command, env=env)
