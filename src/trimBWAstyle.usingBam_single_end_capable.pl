#!/usr/bin/env perl

# AUTHOR: Joseph Fass
# LAST REVISED: January 2010
# Hacked: Rogan Carr rogan@uw.edu
# Date: November 22, 2013
#  Added single-end capability
# Revised: April 30th, 2014
#  Fixed
#   The off-by-1 error in the last-good-base
#   Returning a length of 1 for a bad read
#   Trimming based on cumulative scores
# Hacked: Alex Eng engal@uw.edu
# Date: November 27, 2017
#  Changed automatic output suffixes to be individually specified file names
# The Bioinformatics Core at UC Davis Genome Center
# http://bioinformatics.ucdavis.edu
# Copyright (c) 2010 The Regents of University of California, Davis Campus.
# All rights reserved.

use strict;
use Getopt::Std;

my $usage = "\nusage: $0 -f <unmapped bam file> [-q # -o # -l #]\n\n".
  "Trims sequences from BAM file a la Heng Li's bwa '-q' option\n\n".
  "\t-f    \tUnmapped bam file to trim\n".
  "\t-q int\tQuality threshold (default:  20).  To trim q2, then set -q to 3\n".
  "\t-o int\tFastq offset value (default:  33).  Sanger fastq is default, XXchangeXtoX65XX use 33 for Illumina.\n".
  "\t-l int\tLength cutoff of trimmed reads (default:  60).\n".
  "\t      \tIf both pairs are > length cutoff, then they will be put into read1/read2 fastq.\n".
  "\t      \tIf only 1 read is > length cutoff, then it will be put into a singleton file.\n".
  "\t-s Single-End Read [HACK]\n\n".
  "For specifics about quality trimming bwa style (see bwa's bwa_trim_read() function in bwaseqio.c)\n\n";

our($opt_q,$opt_o,$opt_f,$opt_l,$opt_h,$opt_s,$opt_a,$opt_b,$opt_c);
getopts('shq:o:f:l:a:b:c:');
if (!defined($opt_q) or !($opt_q =~ m/^[0-9]+$/)) {$opt_q = 20;}
if (!defined($opt_o) or !($opt_o =~ m/^[0-9]+$/)) {$opt_o = 33;}
if (!defined($opt_l) or !($opt_l =~ m/^[0-9]+$/)) {$opt_l = 60;}
# The output prefix
# if (!defined($opt_p) or ($opt_p eq '')) {
#   $opt_p = $opt_f; # default
#   $opt_p =~ s/.bam$//;
# }

if (!$opt_f || $opt_h) {die $usage;}

chomp(my $samtools = `which samtools`);
if (!$samtools) {
    die "Samtools not found.  Please check path.\n";
}

#my $pos;  my $maxPos;  my $area;  my $maxArea;

my $read1_file = '';
my $read2_file = '';
my $singletonFile = '';

if (defined($opt_a)){
  $read1_file = $opt_a;
}
if (defined($opt_b)){
  $read2_file = $opt_b;
}
if (defined($opt_c)){
  $singletonFile = $opt_c;
}
# $read1_file .= '.trimmed.1.fastq';
# $read2_file .= '.trimmed.2.fastq';
# $singletonFile .= '.trimmed.singleton.fastq';
#$read1_file =~ s/bam/trimmed.1.fastq/;
#$read2_file =~ s/bam/trimmed.2.fastq/;
#$singletonFile =~ s/bam/trimmed.singleton.fastq/;
my $quality_reads = 0;
my $quality_bases = 0;

open (FILE, "$samtools view $opt_f |") or die "Can't open bam for reading, $opt_f:  $!\n";
# HACK for single-end
unless ( $opt_s ) {
  open (READ1, ">$read1_file") or die "Can't open read1 file, $read1_file\n";
  open (READ2, ">$read2_file") or die "Can't open read2 file, $read2_file\n";
}
open (SINGLETON, ">$singletonFile") or die "Can't open singleton file, $singletonFile:  $!\n";

while (<FILE>) {
    chomp;
    my $read1 = $_;
    my @array = split/\t/, $read1;
    my $read1_name = $array[0];
    if ( ! $opt_s ) {
      $read1_name .= "/1";
    }
    my $read1_s = $array[9];
    my $read1_q = $array[10];
    my $read1_length = &checkPos($read1_q);
    my $read1_seq = substr($read1_s,0,$read1_length);  
    my $read1_qual = substr($read1_q,0,$read1_length);

    my ($read2, $read2_name, $read2_s, $read2_q, $read2_length, $read2_seq, $read2_qual);
    # HACK for single-end
    unless ( $opt_s ) {
      $read2 = <FILE>;
      undef(@array);
      @array = split/\t/, $read2;
      $read2_name = $array[0]."/2";
      $read2_s = $array[9];
      $read2_q = $array[10];
      $read2_length = &checkPos($read2_q);
      $read2_seq = substr($read2_s,0,$read2_length);  
      $read2_qual = substr($read2_q,0,$read2_length);
    } else {
      $read2_length = 0;# HACK for single-end
    }
    if ($read1_length >= $opt_l && $read2_length >= $opt_l) {
	print READ1 "\@$read1_name\n$read1_seq\n\+$read1_name\n$read1_qual\n";
	print READ2 "\@$read2_name\n$read2_seq\n\+$read2_name\n$read2_qual\n";
	$quality_reads += 2;
	$quality_bases += $read1_length;
	$quality_bases += $read2_length;
    } elsif ($read1_length < $opt_l && $read2_length >= $opt_l) {
	print SINGLETON "\@$read2_name\n$read2_seq\n\+$read2_name\n$read2_qual\n";
	$quality_reads++;
	$quality_bases += $read2_length;
    } elsif ($read2_length < $opt_l && $read1_length >= $opt_l) { # Single-end always takes this branch
	print SINGLETON "\@$read1_name\n$read1_seq\n\+$read1_name\n$read1_qual\n";
	$quality_reads++;
	$quality_bases += $read1_length;
    }

}

print "Quality reads: $quality_reads\n";
print "Quality bases: $quality_bases\n";

close FILE;
unless ( $opt_s ) {
  close READ1;
  close READ2;
}
close SINGLETON;

print STDERR "Done\n";

# sub checkPos {
#     my $q = shift;
#     $pos = length($q);
#     $maxPos = $pos;
#     $area = 0;
#     $maxArea = 0;
#     while ($pos>0 && $area>=0) { # This is cumulative, so it's not necessarily exiting on the first good base
# 	$area += $opt_q - (ord(substr($q,$pos-1,1))-$opt_o);
# 	if ($area > $maxArea) {
# 	    $maxArea = $area;
# 	    $maxPos = $pos; # Off by 1 error: returns the position of the last bad base
# 	}
# 	$pos--;
	
#     }  
#     $maxPos = ($pos == 0 ? 1 : $maxPos);
#     return $maxPos;
# }

sub checkPos {
    my $q = shift;
    my $pos = length($q);
    my $length = $pos;
    my $score=-1;
    
    # Count from the end of the read towards the beginning
    #  Exit if we hit the beginning, or if we get a good score.
    while ( $pos>0 && $score < 0 ) {
      # Get the quality for this base
      my $quality = ord(substr($q,$pos-1,1));
      
      # Calculate the score above the threshold
      $score = ($quality-$opt_o) - $opt_q;
      
      if ($score < 0) {
	$length = $pos-1;
      }
      $pos--;
    }
    $length = ($pos == 0 ? 0 : $length);
    return $length;
}

