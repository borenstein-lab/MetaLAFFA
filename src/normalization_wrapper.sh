#!/bin/bash
#
# Author: Alex Eng
# Date: 11/29/2017
#
# Bash wrapper script to normalize KO profiles
#
# Usage:
# normalization_wrapper.sh [-h] [--musicc musicc] [--musicc_correct musicc_correct] ko_profiles normalization_method output
#
# Arguments:
# ko_profiles						: KO profile table to normalize
# normalization_method				: Method to use to normalize the KO profiles (none, musicc)
# output 							: Output file for normalized KO profiles
#
# Options:
# -h								: Print this help information and exit
# --musicc_env musicc_env					: Path to MUSiCC virtual environment (default: /net/borenstein/vol1/PROGRAMS/virtualenv-15.1.0/musicc_env/)
# --musicc_correct musicc_correct	: MUSiCC abundance correction method to use (use_generic, learn_model)

# Setting constants
NORMALIZATION_METHODS=(none musicc)
MUSICC_CORRECTION_METHODS=(use_generic learn_model)

# Initialize variables that need to be set for normalization
ko_profiles=""
normalization_method=""
output=""
musicc_env=/net/borenstein/vol1/PROGRAMS/virtualenv-15.1.0/musicc_env/
musicc_correct=""

# Parse arguments
position_args=()
while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		-h)
			printf "%s\n\n" "Bash wrapper script to normalize KO profiles"
			printf "%s\n" "Usage:"
			printf "%s\n\n" "normalization_wrapper.sh [-h] [--musicc musicc] [--musicc_correct musicc_correct] ko_profiles normalization_method output"
			printf "%s\n" "Arguments:"
			printf "%-32s%s\n" "ko_profiles" ": KO profile table to normalize"
			printf "%-32s%s\n" "normalization_method" ": Method to use to normalize the KO profiles (none, musicc)"
			printf "%-32s%s\n" "output" ": Output file for normalized KO profiles"
			printf "\n"
			printf "%s\n" "Options:"
			printf "%-32s%s\n" "-h" ": Print this help information and exit"
			printf "%-32s%s\n" "--musicc_env musicc_env" ": Path to MUSiCC virtual environment (default: /net/borenstein/vol1/PROGRAMS/virtualenv-15.1.0/musicc_env/)"
			printf "%-32s%s\n" "--msuicc_correct msuicc_correct" ": MUSiCC abundance correction method to use (use_generic, learn_model)"
			exit
			;;
		--musicc_env)
			if [ -e $2 ]
			then
				musicc_env=$2
				shift
			else
				(>&2 echo "The specified file for musicc does not exist (${2})")
				exit 1
			fi
			;;
		--musicc_correct) musicc_correct=$2; shift;;
		*) position_args+=($1);;
	esac
shift
done

# Check if there are enough required positional arguments
if [ ${#position_args[@]} -lt 3 ]
then
	(>&2 echo "Missing one or more required arguments")
	exit 1
fi

# Set positional argument variables
ko_profiles=${position_args[0]}
normalization_method=${position_args[1]}
output=${position_args[2]}

# Check if the required ko profiles file exists
if [ ! -e $ko_profiles ]
then
	(>&2 echo "The specified ko profiles file does not exist (${ko_profiles})")
	exit 1
fi

# Check if the normalization method is recognized
normalization_method_recognized=false
for method in "${NORMALIZATION_METHODS[@]}"
do
	if [ $normalization_method == $method ]
	then
		normalization_method_recognized=true
	fi
done
if [ $normalization_method_recognized == false ]
then
	(>&2 echo "The specified normalization method is not recognized (${normalization_method})")
	exit 1
fi

# Perform normalization based on the method chosen
case $normalization_method in

	# If none, we don't perform any normalization and just write the input to the output
	none) zcat $ko_profiles > $output;;

	# If musicc, then we run MUSiCC
	musicc)
		
		# Check if the abundance correction method is recognized
		musicc_correct_recognized=false
		if [ ! -z $musicc_correct ]
		then
			for method in "${MUSICC_CORRECTION_METHODS[@]}"
			do
				if [ $musicc_correct == $method ]
				then
					musicc_correct_recognized=true
				fi
			done
			if [ $musicc_correct_recognized == false ]
			then
				(>&2 echo "The specified MUSiCC abundance correction method is not recognized (${musicc_correct})")
				exit 1
			fi
		else
			(>&2 echo "If running MUSiCC, an abundance correction method must be specified")
			exit 1
		fi

		# Run MUSiCC
		source ${musicc_env}/bin/activate
		run_musicc.py -n -c $musicc_correct -o $output $ko_profiles
		;;
esac