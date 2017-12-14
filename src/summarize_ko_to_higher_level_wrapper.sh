#!/bin/bash
#
# Author: Alex Eng
# Date: 11/29/2017
#
# Bash wrapper script to summarize KO profiles to higher levels (e.g. module, pathway, etc.)
#
# Usage:
# summarize_ko_to_higher_level_wrapper.sh [--simple_mapper --simple_mapper] [-h] ko_profiles summary_method mapping_matrix output
#
# Arguments:
# ko_profiles                   : KO profile table to summarize to a higher functional level
# summary_method                : Method to use to summarize the KO profiles (fractional, whole)
# mapping_matrix                : Matrix mapping KOs to higher functional summary level
# functional_level				: Name of the functional level being summarized to
# output                        : Output file for higher-level functional profiles
#
# Options:
# --simple_mapper simple_mapper	: Location of simple mapper program (default: src/summarize_ko_to_higher_level.py)
# -h                            : Print this help information and exit

# Setting constants
SUMMARY_METHODS=(fractional whole)

# Initialize variables that need to be set for summarizing
ko_profiles=""
summary_method=""
mapping_matrix=""
functional_level=""
output=""
simple_mapper=src/summarize_ko_to_higher_level.py

# Parse arguments
position_args=()
while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		--simple_mapper)
			if [ -e $2 ]
			then
				simple_mapper=$2
				shift
			else
				(>&2 echo "The specified file for simple_mapper does not exist (${2})")
				exit 1
			fi
			;;
		-h)
			printf "%s\n\n" "Bash wrapper script to summarize KO profiles to higher levels (e.g. module, pathway, etc.)"
			printf "%s\n" "Usage:"
			printf "%s\n\n" "summarize_ko_to_higher_level_wrapper.sh [--simple_mapper --simple_mapper] [-h] ko_profiles summary_method mapping_matrix output"
			printf "%s\n" "Arguments:"
			printf "%-30s%s\n" "ko_profiles" ": KO profile table to summarize to a higher functional level"
			printf "%-30s%s\n" "summary_method" ": Method to use to summarize the KO profiles (fractional, whole)"
			printf "%-30s%s\n" "mapping_matrix" ": Matrix mapping KOs to higher functional summary level"
			printf "%-30s%s\n" "functional_level" ": Name of the functional level being summarized to"
			printf "%-30s%s\n" "output" ": Output file for higher-level functional profiles"
			printf "\n"
			printf "%s\n" "Options:"
			printf "%-30s%s\n" "--simple_mapper" ": Location of simple mapper program (default: src/summarize_ko_to_higher_level.py)"
			printf "%-30s%s\n" "-h" ": Print this help information and exit"
			exit
			;;
		*) position_args+=($1);;
	esac
shift
done

# Check if there are enough required positional arguments
if [ ${#position_args[@]} -lt 4 ]
then
	(>&2 echo "Missing one or more required arguments")
	exit 1
fi

# Set positional argument variables
ko_profiles=${position_args[0]}
summary_method=${position_args[1]}
mapping_matrix=${position_args[2]}
functional_level=${position_args[3]}
output=${position_args[4]}

# Check if the required ko profiles file exists
if [ ! -e $ko_profiles ]
then
	(>&2 echo "The specified ko profiles file does not exist (${ko_profiles})")
	exit 1
fi

# Check if the summary method is recognized
summary_method_recognized=false
for method in "${SUMMARY_METHODS[@]}"
do
	if [ $summary_method == $method ]
	then
		summary_method_recognized=true
	fi
done
if [ $summary_method_recognized == false ]
then
	(>&2 echo "The specified summary method is not recognized (${summary_method})")
	exit 1
fi

# Check if the mapping matrix exists
if [ ! -e $mapping_matrix ]
then
    (>&2 echo "The specified mapping matrix is not recognized (${mapping_matrix})")
    exit 1
fi

# Perform functional summary based on the method chosen
case $summary_method in

	# If fractional, we use the simple mapping program with the fractional mapping setting
	fractional) $simple_mapper $ko_profiles $summary_method $mapping_matrix $functional_level > $output;;

    # If whole, we use the simple mapping program with the whole mapping setting
    whole) $simple_mapper $ko_profiles $summary_method $mapping_matrix $functional_level > $output;;
esac
