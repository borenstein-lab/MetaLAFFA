import argparse
from config import env
import subprocess
import config.operation as op

parser = argparse.ArgumentParser(description="Wrapper script to run the pipeline, ensuring that the shell environment for Snakemake is correct.")
parser.add_argument("--submission_wrapper", "-s", default="python3 src/sge_submission_wrapper.py", help="The wrapper script to use to parse job settings and resource requests and submit a job to the cluster (default %(default)s)")
parser.add_argument("--jobscript", "-j", default="src/configured_jobscript.sh", help="The job script to use when an individual cluster job (default %(default)s)")
parser.add_argument("--number_of_jobs", "-n", type=int, default=50, help="Number of jobs have submitted to the cluster queue at any one time, or if running locally, the number of cores to use for running jobs in parallel (default %(default)s).")
parser.add_argument("--wait", "-w", type=int, default=60, help="Number of seconds to wait after a job finishes before checking that the output exists, this can avoid Snakemake incorrectly marking a step as failed when a file might not be present due to latency over a shared file system (default %(default)s).")
parser.add_argument("--dryrun", "-d", action="store_true", help="If used, determine what steps would be performed and report to the user, but do not actually run the steps.")
parser.add_argument("--local", "-l", action="store_true", help="If used, run metaLAFFA locally rather than on a cluster.")

args = parser.parse_args()

# Create the snakemake command and add options to it as specified by the arguments
command = [op.snakemake]
if not args.local:
    command += ["-c", args.submission_wrapper]
if not args.local:
    command += ["--js", args.jobscript]
command += ["-j", str(args.number_of_jobs)]
command += ["--latency-wait", str(args.wait)]
if args.dryrun:
    command += ["-n"]
command += ["--verbose", "-p"]

# Run the snakemake command
subprocess.run(command, env=env)
