import argparse
import subprocess
import re
import os
import config.operation as op

parser = argparse.ArgumentParser(description="Creates a new metaLAFFA project directory that utilizes the same resources as your original metaLAFFA installation, but with its own configuration files that you can customize to the new project.")
parser.add_argument("new_directory", help="Path to new project directory.")

args = parser.parse_args()

# Copy the base configuration files and the pipeline wrapper
subprocess.run(["cp", "-r", "config", args.new_directory])
subprocess.run(["cp", op.pipeline_step_list, args.new_directory])
subprocess.run(["cp", "metaLAFFA.py", args.new_directory])


def modify_config_module(original_module_name, paths_to_modify):
    """
    Modifies relative paths so they are absolute paths by prefixing them with the path to this directory.

    :param original_module_name: Name of original configuration module to modify
    :param paths_to_modify: List of path variables to modify
    :return: None
    """
    new_module_name = original_module_name + ".tmp"

    with open(original_module_name) as original_module, open(new_module_name) as new_module:
        for line in original_module:
            fields = line.strip().split(" ")
            if fields[0] in paths_to_modify:
                if re.match("^/", fields[2]) is None:
                    fields[2] = os.getcwd() + "/" + fields[2]
            new_module.write(" ".join(fields) + "\n")

    subprocess.run(["mv", new_module_name, original_module_name])


# Modify important non-absolute paths to use the resources from the original installation
all_paths_to_modify = ["bitmask_directory", "database_directory", "index_directory", "gene_normalization_directory", "gene_to_ortholog_directory", "ortholog_to_grouping_directory", "source_directory", "snakefile"]

modify_config_module("config/file_organization.py", all_paths_to_modify)
modify_config_module("config/operation.py", all_paths_to_modify)
