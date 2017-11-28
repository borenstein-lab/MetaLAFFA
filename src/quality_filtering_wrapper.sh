#!/bin/bash
#
# Author: Alex Eng
# Date: 11/27/2017
#
# Bash wrapper script to run trimBWAstyle.usingBam_single_end_capable.pl for trimming and filtering reads by quality from fastq files
#
# Usage:
# quality_filtering_wrapper.sh [--fastq_to_sam fastq_to_sam] [--fix_paired_fastq fix_paired_fastq] [-h] [--paired_fastq paired_fastq] [--paired_fastq_output paired_fastq_output] [--quality_format quality_format] [--singleton_output singleton_output] [--sort_order sort_order] [--trim_bwa_style trim_bwa_style] fastq sample_name output
#
# Arguments:
# fastq										: FASTQ file of reads to filter
# sample_name								: Name of sample the FASTQ is from
# output									: Output file for filtered FASTQ file
#
# Options:
# --fastq_to_sam fastq_to_sam				: Location of fastq_to_sam program (default: /net/gs/vol3/software/modules-sw/picard/1.111/Linux/RHEL6/x86_64/FastqToSam.jar)
# --fix_paired_fastq fix_paired_fastq		: Location of program to fix paired fastq files for conversion to sam format (default: /net/borenstein/vol1/PIPELINE/Updated_Functional_Annotation_Pipeline/src/fix_paired_fastq.py)
# -h										: Print this help information and exit
# --paired_fastq paired_fastq 				: FASTQ file of paired reads
# --paired_fastq_output paired_fastq_output	: Output file for filtered paired FASTQ file
# --quality_format quality_format			: Quality format for SAM file (default: Standard)
# --singleton_output singleton_output		: Output file for newly generated singleton reads
# --sort_order sort_order					: Sort order for SAM file (default: coordinate)
# --trim_bwa_style trim_bwa_style			: Location of trimBWAstyle program (default: /net/borenstein/vol1/PIPELINE/Updated_Functional_Annotation_Pipeline/src/trimBWAstyle.usingBam_single_end_capable.pl)

# Initialize variables that need to be set for bmtagger
fastq=""
sample_name=""
output=""
fastq_to_sam=/net/gs/vol3/software/modules-sw/picard/1.111/Linux/RHEL6/x86_64/FastqToSam.jar
fix_paired_fastq=/net/borenstein/vol1/PIPELINE/Updated_Functional_Annotation_Pipeline/src/fix_paired_fastq.py
paired_fastq=""
paired_fastq_output=""
quality_format=Standard
singleton_output=""
sort_order=coordinate
trim_bwa_style=/net/borenstein/vol1/PIPELINE/Updated_Functional_Annotation_Pipeline/src/trimBWAstyle.usingBam_single_end_capable.pl

# Parse arguments
position_args=()
while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		--fastq_to_sam)
			if [ -e $2 ]
			then
				fastq_to_sam=$2
				shift
			else
				(>&2 echo "The specified file for fastq_to_sam does not exist (${2})")
				exit 1
			fi
			;;
		--fix_paired_fastq)
			if [ -e $2 ]
			then
				fix_paired_fastq=$2
				shift
			else
				(>&2 echo "The specified file for fix_paired_fastq does not exist (${2})")
				exit 1
			fi
			;;
		-h)
			printf "%s\n\n" "Bash wrapper script to run trimBWAstyle.usingBam_single_end_capable.pl for trimming and filtering reads by quality from fastq files"
			printf "%s\n" "Usage:"
			printf "%s\n\n" "quality_filtering_wrapper.sh [--fastq_to_sam fastq_to_sam] [--fix_paired_fastq fix_paired_fastq] [-h] [--paired_fastq paired_fastq] [--paired_fastq_output paired_fastq_output] [--quality_format quality_format] [--singleton_output singleton_output] [--sort_order sort_order] [--trim_bwa_style trim_bwa_style] fastq sample_name output"
			printf "%s\n" "Arguments:"
			printf "%-42s%s\n" "fastq" ": FASTQ file of reads to filter"
			printf "%-42s%s\n" "sample_name" ": Name of sample the FASTQ is from"
			printf "%-42s%s\n" "output" ": Output file for filtered reads"
			printf "\n"
			printf "%s\n" "Options:"
			printf "%-42s%s\n" "--fastq_to_sam fastq_to_sam" ": Location of fastq_to_sam program (default: /net/gs/vol3/software/modules-sw/picard/1.111/Linux/RHEL6/x86_64/FastqToSam.jar)"
			printf "%-42s%s\n" "--fix_paired_fastq fix_paired_fastq" ": Location of program to fix paired fastq files for conversion to sam format (default: /net/borenstein/vol1/PIPELINE/Updated_Functional_Annotation_Pipeline/src/fix_paired_fastq.py)"
			printf "%-42s%s\n" "-h" ": Print this help information and exit"
			printf "%-42s%s\n" "--paired_fastq paired_fastq" ": FASTQ file of paired reads"
			printf "%-42s%s\n" "--paired_fastq_output paired_fastq_output" ": Output file for filtered paired reads"
			printf "%-42s%s\n" "--quality_format" ": Quality format for SAM file (default: Standard)"
			printf "%-42s%s\n" "--singleton_output singleton_output" ": Output file for newly generated singleton reads"
			printf "%-42s%s\n" "--sort_order sort_order" ": Sort order for SAM file (default: coordinate)"
			printf "%-42s%s\n" "--trim_bwa_style trim_bwa_style" ": Location of trimBWAstyle program (default: /net/borenstein/vol1/PIPELINE/Updated_Functional_Annotation_Pipeline/src/trimBWAstyle.usingBam_single_end_capable.pl)"
			exit
			;;
		--paired_fastq)
			if [ -e $2 ]
			then
				paired_fastq=$2
				shift
			else
				(>&2 echo "The specified file for paired_fastq does not exist (${2})")
				exit 1
			fi
			;;
		--paired_fastq_output) paired_fastq_output=$2; shift;;
		--quality_format) quality_format=$2; shift;;
		--singleton_output) singleton_output=$2; shift;;
		--sort_order) sort_order=$2; shift;;
		--trim_bwa_style)
			if [ -e $2 ]
			then
				trim_bwa_style=$2
				shift
			else
				(>&2 echo "The specified file for trim_bwa_style does not exist (${2})")
				exit 1
			fi
			;;
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
fastq=${position_args[0]}
sample_name=${position_args[1]}
output=${position_args[2]}

# Check if the required fastq file exists
if [ ! -e $fastq ]
then
	(>&2 echo "The specified fastq file does not exist (${2})")
	exit 1
fi

# Load the samtools module
module load samtools/latest

# If there is no paired fastq file, run trimBWAstyle in singleton mode
if [ -z $paired_fastq ]
then

	# Convert the fastq to SAM format
	java -jar $fastq_to_sam F1=$fastq O=${fastq}.singleton.sam V=$quality_format SO=$sort_order SM=$sample_name

	# Run trimBWAstyle
	$trim_bwa_style -f ${fastq}.singleton.sam -c $output -s
	rm ${fastq}.singleton.sam

# Otherwise, we run trimBWAstyle in paired mode
else

	# If no paired fastq output file or new singleton output file was specified, exit
	if [ -z $paired_fastq_output ] || [ -z $singleton_output ]
	then 
		echo "If a paired fastq file is specified, a paired fastq output file and a new singleton output file must also be specified"
		exit 1
	fi

	# Fix the paired read files for conversion to SAM format
	$fix_paired_fastq $fastq 1 > ${fastq}.fixed
	$fix_paired_fastq $paired_fastq 2 > ${paired_fastq}.fixed

	# Convert the paired read files to SAM format
	java -jar $fastq_to_sam F1=${fastq}.fixed F2=${paired_fastq}.fixed O=${fastq}.paired.sam V=$quality_format SO=$sort_order SM=$sample_name
	rm ${fastq}.fixed ${paired_fastq}.fixed

	# Run trimBWAstyle
	$trim_bwa_style -f ${fastq}.paired.sam -a $output -b $paired_fastq_output -c $singleton_output
	rm ${fastq}.paired.sam
fi