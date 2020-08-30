#!/usr/bin/env python

import argparse
import subprocess
import config.operation as op
import config.cluster as cl

parser = argparse.ArgumentParser(description="Wrapper script to run the pipeline, ensuring that the shell environment for Snakemake is correct.")
parser.add_argument("--use-cluster", action="store_true", help="MetaLAFFA should use the configured cluster environment for running jobs.")
parser.add_argument("--cluster", "-c", default=cl.submission_wrapper, help="The submission wrapper to use when submitting an individual cluster job (default: %(default)s)")
parser.add_argument("--jobscript", "--js", default=cl.jobscript, help="The job script to use when submitting an individual cluster job (default: %(default)s)")
parser.add_argument("--cores", "--jobs", "-j", type=int, default=50, help="Number of jobs to have submitted to the cluster queue at any one time, or if running locally, the number of cores to use for running jobs in parallel (default: %(default)s).")
parser.add_argument("--latency-wait", "--output-wait", "-w", type=int, default=60, help="Number of seconds to wait after a job finishes before checking that the output exists. This can avoid Snakemake incorrectly marking a step as failed when a file might not be immediately visible due to network latency (default: %(default)s).")

args, remaining_args = parser.parse_known_args()

# Create the snakemake command and add options to it as specified by the arguments
command = [op.snakemake]
if args.use_cluster:
    command += ["-c", args.cluster]
    command += ["--js", args.jobscript]
command += ["-j", str(args.cores)]
command += ["-w", str(args.latency_wait)]
command += remaining_args

# Run the snakemake command
subprocess.run(command)
