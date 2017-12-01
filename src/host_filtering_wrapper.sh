#!/bin/bash
#
# Author: Alex Eng
# Date: 11/27/2017
#
# Bash wrapper script to run bmtagger for filtering out host reads from fastq files
#
# Usage:
# host_filtering_wrapper.sh [--blastdb blastdb] [--bmfiles bmfiles] [--bmfilter bmfilter] [--bmtagger bmtagger] [--extract_fa extract_fa] [-h] [--marked_read_remover marked_read_remover] [--paired_fastq paired_fastq] [--paired_fastq_output paired_fastq_output] [--read_name_matcher read_name_matcher] [--srindex srindex] [--srprism srprism] fastq marked_reads output
#
# Arguments:
# fastq										: FASTQ file of reads to filter
# marked_reads								: File for list of marked reads
# output									: Output file for filtered FASTQ file
#
# Options:
# --blastdb blastdb							: Location of the database used by blastn (default: /net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.fna)
# --bmfiles bmfiles							: Location of bitmask used by bmfilter (default: /net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.bitmask)
# --bmfilter bmfilter						: Location of the bmfilter program (default: /net/borenstein/vol1/PROGRAMS/bmtagger/miniconda2/bin/bmfilter)
# --bmtagger bmtagger						: Location of the bmtagger program (default: /net/borenstein/vol1/PROGRAMS/bmtagger/miniconda2/bin/bmtagger.sh)
# --extract_fa extract_fa					: Location of the extract_fa program (default: src/extract_fullseq)
# -h										: Print this help information and exit
# --marked_read_remover marked_read_remover	: Location of marked read removing program (default: src/remove_marked_reads.py)
# --paired_fastq paired_fastq 				: FASTQ file of paired reads
# --paired_fastq_output paired_fastq_output	: Output file for filtered paired FASTQ file
# --read_name_matcher read_name_matcher		: Location of read name matching program (default: src/match_read_names.py)
# --srindex srindex							: Location of index used by srprism (default: /net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.srprism)
# --srprism srprism							: Location of srprism program (default: src/srprism)

# Initialize variables that need to be set for bmtagger
fastq=""
marked_reads=""
output=""
blastdb=/net/borenstein/vol1/DATA_BLASTDBs/img_core_v400.human_reference.fna
bmfiles=/net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.bitmask
bmfilter=/net/borenstein/vol1/PROGRAMS/bmtagger/miniconda2/bin/bmfilter
bmtagger=/net/borenstein/vol1/PROGRAMS/bmtagger/miniconda2/bin/bmtagger.sh
extract_fa=src/extract_fullseq
marked_read_remover=src/remove_marked_reads.py
paired_fastq=""
paired_fastq_output=""
read_name_matcher=src/match_read_names.py
srindex=/net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.srprism
srprism=src/srprism

# Parse arguments
position_args=()
while [[ $# > 0 ]]
do
	key="$1"
	case $key in
		--blastdb)
			if [ -e $2 ]
			then
				blastdb=$2
				shift
			else
				(>&2 echo "The specified file for blastdb does not exist (${2})")
				exit 1
			fi
			;;
		--bmfiles)
			if [ -e $2 ]
			then
				bmfiles=$2
				shift
			else
				(>&2 echo "The specified file for bmfiles does not exist (${2})")
				exit 1
			fi
			;;
		--bmfilter)
			if [ -e $2 ]
			then
				bmfilter=$2
				shift
			else
				(>&2 echo "The specified file for bmfilter does not exist (${2})")
				exit 1
			fi
			;;
		--bmtagger)
			if [ -e $2 ]
			then
				bmtagger=$2
				shift
			else
				(>&2 echo "The specified file for bmtagger does not exist (${2})")
				exit 1
			fi
			;;
		--extract_fa)
			if [ -e $2 ]
			then
				extract_fa=$2
				shift
			else
				(>&2 echo "The specified file for extract_fa does not exist (${2})")
				exit 1
			fi
			;;
		-h)
			printf "%s\n\n" "Bash wrapper script to run bmtagger for filtering out host reads from fastq files"
			printf "%s\n" "Usage:"
			printf "%s\n\n" "host_filtering_wrapper.sh [--bmfiles bmfiles] [--bmfilter bmfilter] ]--bmtagger bmtagger] [--extract_fa extract_fa] [-h] [--marked_read_remover marked_read_remover] [--paired_fastq paired_fastq] [--paired_fastq_output paired_fastq_output] [--read_name_matcher read_name_matcher] [--srindex srindex] [--srprism srprism] fastq marked_reads output"
			printf "%s\n" "Arguments:"
			printf "%-42s%s\n" "fastq" ": FASTQ file of reads to filter"
			printf "%-42s%s\n" "marked_reads" ": File for list of marked reads"
			printf "%-42s%s\n" "output" ": Output file for filtered reads"
			printf "\n"
			printf "%s\n" "Options:"
			printf "%-42s%s\n" "--blastdb blastdb" ": Location of the database used by blastn (default: /net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.fna)"
			printf "%-42s%s\n" "--bmfiles bmfiles" ": Location of bitmask used by bmfilter (default: /net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.bitmask)"
			printf "%-42s%s\n" "--bmfilter bmfilter" ": Location of the bmfilter program (default: /net/borenstein/vol1/PROGRAMS/bmtagger/miniconda2/bin/bmfilter)"
			printf "%-42s%s\n" "--bmtagger bmtagger" ": Location of the bmtagger program (default: /net/borenstein/vol1/PROGRAMS/bmtagger/miniconda2/bin/bmtagger.sh)"
			printf "%-42s%s\n" "--extract_fa extract_fa" ": Location of the extract_fa program (default: src/extract_fullseq)"
			printf "%-42s%s\n" "-h" ": Print this help information and exit"
			printf "%-42s%s\n" "--marked_read_remover marked_read_remover" ": Location of marked read removing program (default: src/remove_marked_reads.py)"
			printf "%-42s%s\n" "--paired_fastq paired_fastq" ": FASTQ file of paired reads"
			printf "%-42s%s\n" "--paired_fastq_output paired_fastq_output" ": Output file for filtered paired reads"
			printf "%-42s%s\n" "--read_name_matcher read_name_matcher" ": Location of read name matching program (default: src/match_read_names.py)"
			printf "%-42s%s\n" "--srindex srindex" ": Location of index used by srprism (default: /net/borenstein/vol1/PIPELINE/ReadFiltering/HumanReferenceFiles/img_core_v400.human_reference.srprism)"
			printf "%-42s%s\n" "--srprism srprism" ": Location of srprism program (default: src/srprism)"
			exit
			;;
		--marked_read_remover)
			if [ -e $2 ]
			then
				marked_read_remover=$2
				shift
			else
				(>&2 echo "The specified file for marked_read_remover does not exist (${2})")
				exit 1
			fi
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
		--read_name_matcher)
			if [ -e $2 ]
			then
				read_name_matcher=$2
				shift
			else
				(>&2 echo "The specified file for read_name_matcher does not exist (${2})")
				exit 1
			fi
			;;
		--srindex)
			if [ -e $2 ]
			then
				srindex=$2
				shift
			else
				(>&2 echo "The specified file for srindex does not exist (${2})")
				exit 1
			fi
			;;
		--srprism)
			if [ -e $2 ]
			then
				srprism=$2
				shift
			else
				(>&2 echo "The specified file for srprism does not exist (${2})")
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
marked_reads=${position_args[1]}
output=${position_args[2]}

# Check if the required fastq file exists
if [ ! -e $fastq ]
then
	(>&2 echo "The specified fastq file does not exist (${2})")
	exit 1
fi

# If zipped, generate unzipped versions of the required input fastq file
fastq_gzipped=false
if [ $(echo $fastq | grep "\.gz$" | wc -l) -eq 1 ]
then
	fastq_gzipped=true
	zcat $fastq > ${fastq}.unzipped
	fastq=${fastq}.unzipped
fi

# Export environment variables for bmtagger
export BMFILTER=$bmfilter EXTRACT_FA=$extract_fa SRPRISM=$srprism

# If there is no paired end fastq, run bmtagger in singleton mode
if [  -z $paired_fastq ]
then
	$bmtagger -q1 -1 $fastq -o $marked_reads -b $bmfiles -d $blastdb -x $srindex

# Otherwise, make sure read names match properly between files and then run bmtagger in paired read mode
else

	# If no paired fastq output file was specified, exit
	if [ -z $paired_fastq_output ]
	then 
		echo "If a paired fastq file is specified, a paired fastq output file must also be specified"
		exit 1
	else
		$read_name_matcher $fastq $paired_fastq > ${paired_fastq}.matched
		$bmtagger -q1 -1 $fastq -2 ${paired_fastq}.matched -o $marked_reads -b $bmfiles -d $blastdb -x $srindex
		rm ${paired_fastq}.matched
	fi
fi

# If we had to unzip the input fastq, remove the unzipped version
if [ $fastq_gzipped == "true" ]
then
	fastq=$(echo $fastq | sed 's/\(.*\)\.unzipped/\1/')
	rm ${fastq}.unzipped
fi

# Remove the marked reads from the original fastq
$marked_read_remover $marked_reads $fastq > $output

# If there is a paired end fastq, remove the marked reads from the paired file too
if [ ! -z $paired_fastq ]
then
	$marked_read_remover $marked_reads $paired_fastq > $paired_fastq_output
fi