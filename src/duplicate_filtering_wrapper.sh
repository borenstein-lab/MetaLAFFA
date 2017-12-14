#!/bin/bash
#
# Author: Alex Eng
# Date: 11/27/2017
#
# Bash wrapper script to run MarkDuplicates for filtering out duplicate reads from fastq files
#
# Usage:
# duplicate_filtering_wrapper.sh [--duplicate_marker duplicate_marker] [--extract_duplicates extract_duplicates] [--fastq_to_sam fastq_to_sam] [--fix_paired_fastq fix_paired_fastq] [-h] [--java java] [--marked_read_remover marked_read_remover] [--paired_fastq paired_fastq] [--paired_fastq_output paired_fastq_output] [--quality_format quality_format] [--samtools samtools] [--sort_order sort_order] fastq sample_name metric_file marked_reads output
#
# Arguments:
# fastq                                     : FASTQ file of reads to filter
# sample_name                               : Name of sample the FASTQ is from
# metric_file                               : Output file for library complexity metrics
# marked_reads                              : File for list of marked reads
# output                                    : Output file for filtered FASTQ file
#
# Options:
# -- duplicate_marker duplicate_marker      : Location of the program to mark duplicates in a SAM file (default: /net/borenstein/vol1/PROGRAMS/DPWG_ProcessingTools/EstimateLibraryComplexity.jar)
# --extract_duplicates extract_duplicates   : Location of duplicate extraction program (default: src/extract_duplicates.py)
# --fastq_to_sam fastq_to_sam               : Location of fastq_to_sam program (default: /net/gs/vol3/software/modules-sw/picard/1.111/Linux/RHEL6/x86_64/FastqToSam.jar)
# --fix_paired_fastq fix_paired_fastq       : Location of program to fix paired fastq files for conversion to sam format (default: src/fix_paired_fastq.py)
# -h                                        : Print this help information and exit
# --java java                               : Location of the Java binary (default: /net/borenstein/vol1/PROGRAMS/java/jre1.8.0_151/bin/java)
# --marked_read_remover marked_read_remover : Location of marked read removing program (default: src/remove_marked_reads.py)
# --paired_fastq paired_fastq               : FASTQ file of paired reads
# --paired_fastq_output paired_fastq_output : Output file for filtered paired FASTQ file
# --quality_format quality_format           : Quality format for SAM file (default: Standard)
# --samtools samtools                       : Location of samtools program (default: /net/gs/vol3/software/modules-sw/samtools/0.1.8/Linux/RHEL6/x86_64/bin/samtools)
# --sort_order sort_order                   : Sort order for SAM file (default: coordinate)

# Initialize variables that need to be set for bmtagger
fastq=""
sample_name=""
metric_file=""
marked_reads=""
output=""
duplicate_marker=/net/borenstein/vol1/PROGRAMS/DPWG_ProcessingTools/EstimateLibraryComplexity.jar
extract_duplicates=src/extract_duplicates.py
fastq_to_sam=/net/gs/vol3/software/modules-sw/picard/1.111/Linux/RHEL6/x86_64/FastqToSam.jar
fix_paired_fastq=src/fix_paired_fastq.py
java=/net/borenstein/vol1/PROGRAMS/java/jre1.8.0_151/bin/java
marked_read_remover=src/remove_marked_reads.py
paired_fastq=""
paired_fastq_output=""
quality_format=Standard
samtools=/net/gs/vol3/software/modules-sw/samtools/0.1.8/Linux/RHEL6/x86_64/bin/samtools
sort_order=coordinate

# Parse arguments
position_args=()
while [[ $# > 0 ]]
do
    key="$1"
    case $key in
        --duplicate_marker)
            if [ -e $2 ]
            then
                duplicate_marker=$2
                shift
            else
                (>&2 echo "The specified file for duplicate_marker does not exist (${2})")
                exit 1
            fi
            ;;
        --extract_duplicates)
            if [ -e $2 ]
            then
                extract_duplicates=$2
                shift
            else
                (>&2 echo "The specified file for extract_duplicates does not exist (${2})")
                exit 1
            fi
            ;;
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
            printf "%s\n\n" "Bash wrapper script to run MarkDuplicates for filtering out duplicate reads from fastq files"
            printf "%s\n" "Usage:"
            printf "%s\n\n" "duplicate_filtering_wrapper.sh [--extract_duplicates extract_duplicates] [--fastq_to_sam fastq_to_sam] [--fix_paired_fastq fix_paired_fastq] [-h] [--marked_read_remover marked_read_remover] [--paired_fastq paired_fastq] [--paired_fastq_output paired_fastq_output] [--quality_format quality_format] [--samtools samtools] [--sort_order sort_order] fastq sample_name metric_file marked_reads output"
            printf "%s\n" "Arguments:"
            printf "%-42s%s\n" "fastq" ": FASTQ file of reads to filter"
            printf "%-42s%s\n" "sample_name" ": Name of sample the FASTQ is from"
            printf "%-58s%s\n" "metric_file" ": Output file for library complexity metrics"
            printf "%-42s%s\n" "marked_reads" ": File for list of marked reads"
            printf "%-42s%s\n" "output" ": Output file for filtered reads"
            printf "\n"
            printf "%s\n" "Options:"
            printf "%-42s%s\n" "--duplicate_marker duplicate_marker" ": Location of the program to mark duplicates in a SAM file (default: /net/borenstein/vol1/PROGRAMS/DPWG_ProcessingTools/EstimateLibraryComplexity.jar)"
            printf "%-42s%s\n" "--extract_duplicates extract_duplicates" ": Location of duplicate extraction program (default: src/extract_duplicates.py)"
            printf "%-42s%s\n" "--fastq_to_sam fastq_to_sam" ": Location of fastq_to_sam program (default: /net/gs/vol3/software/modules-sw/picard/1.111/Linux/RHEL6/x86_64/FastqToSam.jar)"
            printf "%-42s%s\n" "--fix_paired_fastq fix_paired_fastq" ": Location of program to fix paired fastq files for conversion to sam format (default: src/fix_paired_fastq.py)"
            printf "%-42s%s\n" "-h" ": Print this help information and exit"
            printf "%-42s%s\n" "--java java" ": Location of the Java binary (default: /net/borenstein/vol1/PROGRAMS/java/jre1.8.0_151/bin/java)"
            printf "%-42s%s\n" "--marked_read_remover marked_read_remover" ": Location of marked read removing program (default: src/remove_marked_reads.py)"
            printf "%-42s%s\n" "--paired_fastq paired_fastq" ": FASTQ file of paired reads"
            printf "%-42s%s\n" "--paired_fastq_output paired_fastq_output" ": Output file for filtered paired reads"
            printf "%-42s%s\n" "--quality_format" ": Quality format for SAM file (default: Standard)"
            printf "%-42s%s\n" "--samtools samtools" ": Location of samtools program (default: /net/gs/vol3/software/modules-sw/samtools/0.1.8/Linux/RHEL6/x86_64/bin/samtools)"
            printf "%-42s%s\n" "--sort_order sort_order" ": Sort order for SAM file (default: coordinate)"
            exit
            ;;
        --java)
            if [ -e $2 ]
            then
                java=$2
                shift
            else
                (>&2 echo "The specified file for java does not exist (${2})")
                exit 1
            fi
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
        --quality_format) quality_format=$2; shift;;
        --samtools)
            if [ -e $2 ]
            then
                samtools=$2
                shift
            else
                (>&2 echo "The specified file for samtools does not exist (${2})")
                exit 1
            fi
            ;;
        --sort_order) sort_order=$2; shift;;
        *) position_args+=($1);;
    esac
shift
done

# Check if there are enough required positional arguments
if [ ${#position_args[@]} -lt 5 ]
then
    (>&2 echo "Missing one or more required arguments")
    exit 1
fi

# Set positional argument variables
fastq=${position_args[0]}
sample_name=${position_args[1]}
metric_file=${position_args[2]}
marked_reads=${position_args[3]}
output=${position_args[4]}

# Check if the required fastq file exists
if [ ! -e $fastq ]
then
    (>&2 echo "The specified fastq file does not exist (${2})")
    exit 1
fi

# If there is no paired fastq file, MarkDuplicates will do nothing, so we just copy the input to the output
if [ -z $paired_fastq ]
then
    zcat $fastq > $output

# Otherwise, we run MarkDuplicates
else

    # If no paired fastq output file was specified, exit
    if [ -z $paired_fastq_output ]
    then 
        echo "If a paired fastq file is specified, a paired fastq output file must also be specified"
        exit 1
    fi

    # Fix the paired read files for conversion to SAM format
    # Convert the paired read files to SAM format
    $java -jar $fastq_to_sam F1=$fastq F2=$paired_fastq O=${fastq}.paired.sam V=$quality_format SO=$sort_order SM=$sample_name

    # Run MarkDuplicates
    $java -jar $duplicate_marker I=${fastq}.paired.sam O=${fastq}.elc_output M=$metric_file
    rm ${fastq}.paired.sam

    # Convert the MarkDuplicates output to a parse-able format
    $samtools view ${fastq}.elc_output > ${fastq}.samview
    rm ${fastq}.elc_output

    # Extract the reads that were marked as duplicates
    $extract_duplicates ${fastq}.samview > $marked_reads
    rm ${fastq}.samview

    # Filter the input paired read files
    $marked_read_remover $marked_reads $fastq > $output
    $marked_read_remover $marked_reads $paired_fastq > $paired_fastq_output
fi