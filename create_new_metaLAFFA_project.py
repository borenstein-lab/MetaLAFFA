import argparse
import subprocess
import re
import os
import config.operation as op
cwd = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description="Creates a new metaLAFFA project directory that utilizes the same resources as your original metaLAFFA installation, but with its own configuration files that you can customize to the new project.")
parser.add_argument("new_directory", help="Path to new project directory.")

args = parser.parse_args()

# Get the absolute path to the new project directory (important if calling script from outside original metaLAFFA install directory)
absolute_new_directory = args.new_directory
if re.match("^\"/", absolute_new_directory) is None:
    absolute_new_directory = os.getcwd() + "/" + absolute_new_directory

# Copy the base configuration files and the pipeline wrapper
os.makedirs(args.new_directory, exist_ok=True)
os.makedirs(args.new_directory + "/src", exist_ok=True)

subprocess.run(["cp", "-r", "config", absolute_new_directory + "/config"], cwd=cwd)
subprocess.run(["cp", op.pipeline_step_list, absolute_new_directory], cwd=cwd)
subprocess.run(["cp", "metaLAFFA.py", absolute_new_directory], cwd=cwd)

# Configure the jobscript for the new project
with open(cwd + "/src/template_jobscript.sh") as template, open(args.new_directory + "/src/configured_jobscript.sh", "w") as configured:

    # Identify the location of the python executable, making it an absolute path if it is a relative path
    python_exec_path = op.python
    if re.search("/", python_exec_path) is not None and re.match("^/") is not None:
        python_exec_path = cwd + "/" + python_exec_path
    for line in template:
        configured.write(line.format(python=python_exec_path, metaLAFFA_directory=absolute_new_directory, PYTHONPATH="{{PYTHONPATH}}"))


def modify_config_module(original_module_name, paths_to_modify):
    """
    Modifies relative paths so they are absolute paths by prefixing them with the path to this directory.

    :param original_module_name: Name of original configuration module to modify
    :param paths_to_modify: List of path variables to modify
    :return: None
    """
    new_module_name = original_module_name + ".tmp"

    with open(original_module_name) as original_module, open(new_module_name, "w") as new_module:
        for line in original_module:
            fields = line.strip().split(" ")
            if fields[0] in paths_to_modify:
                if re.match("^\"/", fields[2]) is None:
                    fields[2] = "\"" + os.getcwd() + "/" + fields[2].strip("\"") + "\""
            new_module.write(" ".join(fields) + "\n")

    subprocess.run(["mv", new_module_name, original_module_name])


# Modify important non-absolute paths to use the resources from the original installation
all_paths_to_modify = ["bitmask_directory", "database_directory", "index_directory", "gene_normalization_directory", "gene_to_ortholog_directory", "ortholog_to_grouping_directory", "source_directory", "snakefile"]

modify_config_module(absolute_new_directory + "/config/file_organization.py", all_paths_to_modify)
modify_config_module(absolute_new_directory + "/config/operation.py", all_paths_to_modify)
