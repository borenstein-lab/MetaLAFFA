import config.operation as op
import config.cluster as cl
import config.file_organization as fo
import config.library_functions as lf
import re
import importlib
import copy
import os


def get_working_name(filename):
    """
    Processes the filename and returns the version the pipeline will used based on config settings (e.g. if we are working with zipped files, make sure expected input and output file names are zipped file names).

    :param filename: Filename to convert into its zipped name
    :return: The zipped version of the input file name
    """

    if op.zipped_files:
        if re.search("\\.gz$", filename):
            return filename
        else:
            return filename + ".gz"
    else:
        return filename


def create_input_generator(input_dic):
    """
    Creates a function that takes in a wildcards object and returns the inputs to a pipeline step

    :param input_dic: Dictionary mapping input names to input file patterns
    :return: Function that processes wildcards to produce inputs
    """

    copied_input_dic = copy.deepcopy(input_dic)

    # Create the input generator function
    def input_generator(wildcards):
        # If there are wildcards, create a new copy of the dictionary to modify and replace wildcards in input file name patterns in the new dictionary
        dic_to_return = copy.deepcopy(copied_input_dic)
        if wildcards is not None:
            for input_name in dic_to_return:

                # If there is a single input pattern, format it with the given wildcards
                if type(dic_to_return[input_name]) == str:
                    dic_to_return[input_name] = dic_to_return[input_name].format(wildcards=wildcards)

                # Otherwise, format each string in the input list
                else:
                    formatted_inputs = []
                    for input_pattern in dic_to_return[input_name]:
                        formatted_inputs.append(input_pattern.format(wildcards=wildcards))
                    dic_to_return[input_name] = formatted_inputs

        return dic_to_return

    return input_generator


def parse_pipeline_steps():
    """
    Parses the pipeline step file to determine which steps should be run, what the inputs are to each step, and which steps produce the final desired outputs of the pipeline.

    :return: Dictionary containing information about pipeline steps to run
    """

    step_info = {"FINAL_OUTPUTS": [], "SUMMARIES": []}
    with open(op.pipeline_step_list) as step_file:

        for line in step_file:

            # Only parse uncommented, non-empty lines
            if not re.match("#", line) and line.strip() != "":

                # Split the line on ":" to get the step name and its inputs
                step_name = line.strip().split(":")[0]
                inputs = line.strip().split(":")[1].split(",")

                # If the step name begins with a "$", then it is a desired final output
                if re.match("\\$", step_name):
                    step_name = step_name[1:]
                    step_info["FINAL_OUTPUTS"].append(step_name)

                # If the step name begins with a "*", then it produces a summary statistic file
                if re.match("\\*", step_name):
                    step_name = step_name[1:]
                    step_info["SUMMARIES"].append(step_name)

                step_info[step_name] = inputs

    return step_info


def assign_step_prefix_components(step_name, step_info, step_params):
    """
    Using information about which pipeline steps are used as inputs for other pipeline steps, assign the provenance-encoding prefix is for the output of a given step in the parameter dictionary, recursively assigning prefixes to prior steps when unknown

    :param step_name: Name of step to get prefix for
    :param step_info: Dictionary of prior steps for each step in the pipeline
    :param step_params: Dictionary of step parameters that contain prefixes from prior steps and that prefixes will be assigned to
    :return: None (modifies step_params in place)
    """

    # If the prefix components have not been set for this step yet, set them, otherwise do nothing
    if "prefix_components" not in step_params[step_name]:

        step_module = importlib.import_module(".".join(["config", "steps", step_name]))
        prior_steps = step_info[step_name]

        # The prefix components begin with the base directory for the kind of output (data output or summary output)
        parent_directory = ""
        # If this step produces summary output, then the parent directory is within the summary output directory
        if step_name in step_info["SUMMARIES"]:
            parent_directory += fo.summary_directory

        # Otherwise, the parent directory within the the output directory
        else:
            parent_directory += fo.output_directory

        # Next add the step's specific output directory
        parent_directory += step_module.step_prefix + "/"
        prefix_components = [parent_directory]

        # Gather all of the prefix components for prior steps so that they can be combined later
        for prior_step in prior_steps:

            # If the prior step is SUMMARIES, gather the prefixes for all steps that produce summary tables
            if prior_step == "SUMMARIES":
                for summary_step in step_info["SUMMARIES"]:

                    # If the summary step has prefix components available, add those (skipping the first, which is the parent directory for that step)
                    if "prefix_components" in step_params[summary_step]:
                        prefix_components += step_params[summary_step]["prefix_components"][1:]

                    # Otherwise recursively determine what the prefix components for the prior step are and then add those (skipping the first, which is the parent directory for that step)
                    else:
                        assign_step_prefix_components(summary_step, step_info, step_params)
                        prefix_components += step_params[summary_step]["prefix_components"][1:]

            # Otherwise, skip INPUT prior steps and add prefix components for everything else
            elif prior_step != "INPUT":

                # If the prior step has prefix components available, add those (skipping the first, which is the parent directory for that step)
                if "prefix_components" in step_params[prior_step]:
                    prefix_components += step_params[prior_step]["prefix_components"][1:]

                # Otherwise recursively determine what the prefix components for the prior step are and then add those (skipping the first, which is the parent directory for that step)
                else:
                    assign_step_prefix_components(prior_step, step_info, step_params)
                    prefix_components += step_params[prior_step]["prefix_components"][1:]

        # Ordering should be maintained in prefix components, while duplicates should be removed
        prefix_component_set = set()
        prefix_component_index = 0

        # Consider each prefix component, removing later duplicates, until there are no more duplicates
        while prefix_component_index < len(prefix_components):

            prefix_component = prefix_components[prefix_component_index]

            # If the current component hasn't bee seen yet, add it to the set of seen components and look at the next component
            if prefix_component not in prefix_component_set:
                prefix_component_set.add(prefix_component)
                prefix_component_index += 1

            # Otherwise, remove the current component from the list
            else:
                del prefix_components[prefix_component_index]

        # Create the full output prefix for this step using the base output prefix and operating parameters and add it to the prefix components
        operating_parameters = []
        for parameter in step_module.operating_params:
            operating_parameters.append(re.sub("\\s", "_", parameter))
            operating_parameters.append(re.sub("\\s", "_", str(step_module.operating_params[parameter])))
        prefix_components.append("_".join([step_module.step_prefix] + operating_parameters))

        # Set the prefix components for the step in the step parameter dictionary
        step_params[step_name]["prefix_components"] = prefix_components


def create_step_params(step_info):
    """
    Uses settings from the config module and the desired set of pipeline steps to run to define the parameters for each pipeline step.

    :param step_info: Dictionary defining the inputs of each pipeline step and the final desired outputs
    :return: Dictionary defining the parameters for each pipeline step that will be run.
    """

    # Initialize dictionary that will hold all of the parameters for each step
    step_params = {"raw_final_outputs": set(), "final_outputs": set()}

    # Get the list of step names (excluding special keys in the step information dictionary
    step_names = []
    for step_name in step_info.keys():
        if step_name not in ["FINAL_OUTPUTS", "SUMMARIES"]:
            step_names.append(step_name)

    # Add each step to the step parameter dictionary
    for step_name in step_names:
        step_params[step_name] = {}

    # Go through desired steps in the pipeline and add parameters to the dictionary
    for step_name in step_names:

        # Import the config submodule for this step
        step_module = importlib.import_module(".".join(["config", "steps", step_name]))

        # Check for custom cluster resource requests and use defaults when not specified
        memory = cl.default_time
        if "memory" in step_module.cluster_params:
            memory = step_module.cluster_params["memory"]
        time = cl.default_time
        if "time" in step_module.cluster_params:
            time = step_module.cluster_params["time"]
        options = cl.default_options
        if "options" in step_module.cluster_params:
            options = step_module.cluster_params["options"]
        cores = cl.default_cores
        if "cores" in step_module.cluster_params:
            cores = step_module.cluster_params["cores"]
            step_params[step_name]["threads"] = step_module.cluster_params["cores"] * op.cpu_to_thread_multiplier
        step_params[step_name]["cluster_params"] = cl.sge_cluster_param_string_generator(memory, time, cores, options)

        # Assign the prefix components to this step
        assign_step_prefix_components(step_name, step_info, step_params)
        prefix_components = step_params[step_name]["prefix_components"]

        # Set the output file name patterns
        output_list = step_module.output_list
        for output_index in range(len(output_list)):
            output_list[output_index] = get_working_name(prefix_components[0] + fo.provenance_separator.join(prefix_components[1:]) + "/" + output_list[output_index])

            # Add the output file to the list of raw final output files if indicated
            if step_name in step_info["FINAL_OUTPUTS"]:
                raw_final_output_files = lf.generate_possible_patterns_from_restricted_wildcards([output_list[output_index]], op.wildcard_restrictions)
                for raw_final_output_file in raw_final_output_files:
                    step_params["raw_final_outputs"].add(raw_final_output_file)
        step_params[step_name]["output"] = output_list

        # Set the benchmark file name pattern
        step_params[step_name]["benchmark"] = get_working_name(fo.benchmark_dir + prefix_components[0] + fo.provenance_separator.join(prefix_components[1:]) + "/" + step_module.benchmark_file)

        # Set the code that runs
        step_params[step_name]["rule_function"] = step_module.rule_function

    # Now that step outputs have been determined, we can set step inputs
    for step_name in step_names:

        # Import the config submodule for this step
        step_module = importlib.import_module(".".join(["config", "steps", step_name]))

        # Set the input file name generator function
        prior_steps = step_info[step_name]
        input_dic = step_module.input_dic
        processed_input_dic = {}
        for input_name in input_dic:

            # Get which prior step this input file is from
            prior_step = prior_steps[input_dic[input_name]["prior_step"]]

            # If the prior step is INPUT, then we use the initial data prefix
            if prior_step == "INPUT":
                processed_input_dic[input_name] = get_working_name(fo.initial_data_directory + input_dic[input_name]["file"])

            # Otherwise, if the prior step is not SUMMARIES, we use the prefix components from the prior step to determine the input file name pattern (we wait for SUMMARIES until we've determined each step's output files)
            elif prior_step == "SUMMARIES":
                processed_input_dic[input_name] = []
                for summary_step_name in step_info["SUMMARIES"]:
                    for summary_file in step_params[summary_step_name]["output"]:
                        processed_input_dic[input_name].append(get_working_name(summary_file))

            # Otherwise, we use the prefix components from the prior step to determine the input file name pattern
            else:
                prefix_components = step_params[prior_step]["prefix_components"]
                processed_input_dic[input_name] = get_working_name(prefix_components[0] + fo.provenance_separator.join(prefix_components[1:]) + "/" + input_dic[input_name]["file"])

        # Add the input generator function to the step parameters
        step_params[step_name]["input"] = create_input_generator(processed_input_dic)

    # After all other processing, we can also now set the final processed output files
    for raw_final_output in step_params["raw_final_outputs"]:
        step_params["final_outputs"].add(lf.process_final_output_name(raw_final_output))
        
        # Set the rule for the final output processing step
        step_params["process_final_output"] = {}
        step_params["process_final_output"]["rule_function"] = importlib.import_module(".".join(["config", "steps", "process_final_output"])).rule_function
    return step_params


def make_directories():
    """
    If output and logging directories specified by the configuration module are not present, create them.

    :return: None
    """

    for required_dir in [fo.output_directory, fo.benchmark_dir, fo.log_directory, fo.summary_directory]:
        if not os.path.isdir(required_dir):
            os.makedirs(required_dir)
