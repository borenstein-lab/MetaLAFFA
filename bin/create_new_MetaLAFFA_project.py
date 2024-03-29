#!/usr/bin/env python

import argparse
import subprocess
import re
import os
import config.file_organization as fo
import config.operation as op
import config.cluster as cl

parser = argparse.ArgumentParser(description="Creates a new MetaLAFFA project directory that utilizes the same reference resources as your original MetaLAFFA installation, but with its own configuration files that you can customize to the new project.")
parser.add_argument("new_directory", help="Path to new project directory.")

args = parser.parse_args()

# Get the absolute path to the new project directory
absolute_new_directory = args.new_directory
if re.match("^\"/", absolute_new_directory) is None:
    absolute_new_directory = os.getcwd() + "/" + absolute_new_directory

# Copy the base configuration files and the pipeline wrapper
os.makedirs(args.new_directory, exist_ok=True)
os.makedirs(args.new_directory + "/src", exist_ok=True)

# Create project-specific data directories
for required_directory in fo.required_project_directories:
    if not os.path.isdir(absolute_new_directory + "/" + required_directory):
        os.makedirs(absolute_new_directory + "/" + required_directory)

subprocess.run(["cp", "-r", fo.installation_directory + "../lib/python3.6/config", absolute_new_directory + "/config"])
subprocess.run(["cp", fo.installation_directory + "/" + op.pipeline_step_list, absolute_new_directory])
subprocess.run(["cp", fo.installation_directory + "/" + op.snakefile, absolute_new_directory])
subprocess.run(["cp", fo.installation_directory + "/MetaLAFFA.py", absolute_new_directory])
subprocess.run(["chmod", "+x", absolute_new_directory + "/MetaLAFFA.py"])
subprocess.run(["cp", fo.installation_directory + "/" + cl.submission_wrapper, absolute_new_directory + "/" + cl.submission_wrapper])
subprocess.run(["chmod", "+x", absolute_new_directory + "/" + cl.submission_wrapper])

# Configure the jobscript for the new project
with open(fo.source_directory + "/template_jobscript.sh") as template, open(args.new_directory + "/" + cl.jobscript, "w") as configured:
    for line in template:
        configured.write(line.format(MetaLAFFA_directory=absolute_new_directory, conda_environment=op.conda_env))
