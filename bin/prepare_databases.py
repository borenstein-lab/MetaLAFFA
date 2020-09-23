#!/usr/bin/env python

import argparse
import sys
import os
import subprocess
import glob
import re
import config.file_organization as fo
import config.operation as op
import config.steps

parser = argparse.ArgumentParser(description="Downloads and formats default databases.")
parser.add_argument("--human_reference",
                    "-hr",
                    action="store_true",
                    help="If used, download the human reference + decoy sequences used in the 1000 genomes project (hs37d5).")
parser.add_argument("--ko_mappings",
                    "-km",
                    action="store_true",
                    help="If used, download bacterial ko-to-module and ko-to-pathway mappings from the 2013 version of the KEGG database.")
parser.add_argument("--uniprot",
                    "-u",
                    action="store_true",
                    help="If used, download the default UniProt gene sequence reference database for read mapping.")
parser.add_argument("--force",
                    "-f",
                    action="store_true",
                    help="Normally, if the expected output files already exists, this script will skip generating them. Use this option to force the script re-download and/or regenerate all specified databases, even if they already exist.")

args = parser.parse_args()

# Create any missing directories
for required_directory in fo.required_reference_directories:
    if not os.path.isdir(required_directory):
        os.makedirs(required_directory)

processing_error = False

processing_results = {}
for database in ["human_reference", "ortholog_mapping", "uniprot"]:
    processing_results[database] = {}
    processing_results[database]["error"] = False
    processing_results[database]["stdout"] = fo.database_directory + database + "_processing.out"
    processing_results[database]["stderr"] = fo.database_directory + database + "_processing.err"

sys.stdout.write("Beginning processing of reference data.\n\n")
# Try to create and process default database files
if args.human_reference:
    source_exists = True
    if not os.path.isfile(op.host_database_file) or args.force:
        sys.stdout.write("Downloading human reference.\n")
        try:
            subprocess.run(["wget",
                            "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/%s.fa.gz" % op.host_database,
                            "-P",
                            fo.database_directory],
                           stdout=open(processing_results["human_reference"]["stdout"], "a"),
                           stderr=open(processing_results["human_reference"]["stderr"], "a"))
            subprocess.run(["gunzip",
                            op.host_database_file + ".gz"],
                           stdout=open(processing_results["human_reference"]["stdout"], "a"),
                           stderr=open(processing_results["human_reference"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["human_reference"]["error"] = True
            sys.stderr.write("Error: There was a problem downloading the human reference database in preparation for processing. Please see %s and %s for processing logs.\n" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"]))
            source_exists = False
    else:
        sys.stdout.write("Human reference exists, skipping download.\n")

    if (len(glob.glob(op.host_index + "*.bt2")) == 0 and source_exists) or (source_exists and args.force):
        sys.stdout.write("Creating human reference Bowtie 2 index.\n")
        try:
            # Tag output index files with a ".tmp" tag and remove it once all index files have been generated
            # just in case something interrupts the indexing process
            subprocess.run([config.op.bowtie2_build,
                            op.host_database_file,
                            op.host_index + ".tmp"],
                           stdout=open(processing_results["human_reference"]["stdout"], "a"),
                           stderr=open(processing_results["human_reference"]["stderr"], "a"))
            for temp_index_file in glob.glob(op.host_index + "*.bt2"):
                permanent_file = re.sub("\\.tmp", "", temp_index_file)
                subprocess.run(["mv",
                                temp_index_file,
                                permanent_file],
                               stdout=open(processing_results["human_reference"]["stdout"], "a"),
                               stderr=open(processing_results["human_reference"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["human_reference"]["error"] = True
            sys.stderr.write("Error: There was a problem creating an index for the human reference database. Please see %s and %s for processing logs.\n" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"]))
    else:
        sys.stdout.write("Human reference Bowtie 2 index exists, skipping database indexing.\n")
sys.stdout.write("\n")

if args.ko_mappings:
    missing_mappings = []
    for ortholog_to_grouping_mapping in op.ortholog_to_grouping_mappings:
        if not os.path.isfile(fo.ortholog_to_grouping_directory +
                              ortholog_to_grouping_mapping +
                              op.ortholog_to_grouping_suffix) or args.force:
            missing_mappings.append(ortholog_to_grouping_mapping)
    if len(missing_mappings) > 0:
        sys.stdout.write("Downloading ortholog-to-grouping mappings.\n")
        for ortholog_to_grouping_mapping in missing_mappings:
            if not os.path.isfile(fo.ortholog_to_grouping_directory +
                                  ortholog_to_grouping_mapping +
                                  op.ortholog_to_grouping_suffix) or args.force:
                try:
                    subprocess.run(["wget",
                                    "https://github.com/borenstein-lab/fishtaco/raw/master/fishtaco/data/" + ortholog_to_grouping_mapping,
                                    "-O",
                                    fo.ortholog_to_grouping_directory + ortholog_to_grouping_mapping + op.ortholog_to_grouping_suffix],
                                   stdout=open(processing_results["ortholog_mapping"]["stdout"], "a"),
                                   stderr=open(processing_results["ortholog_mapping"]["stderr"], "a"))
                except (EnvironmentError, subprocess.CalledProcessError):
                    processing_error = True
                    processing_results["ortholog_mapping"]["error"] = True
                    sys.stderr.write("Error: There was a problem downloading the %s ortholog-to-grouping mapping file. Please see %s and %s for processing logs.\n" % (ortholog_to_grouping_mapping, processing_results["ortholog_mapping"]["stdout"], processing_results["ortholog_mapping"]["stderr"]))
    else:
        sys.stdout.write("Ortholog-to-grouping mappings exist, skipping download.\n")
sys.stdout.write("\n")

if args.uniprot:
    source_exists = True
    if not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz") or args.force:
        sys.stdout.write("Downloading UniProt UniRef90 database.\n")
        try:
            subprocess.run(["wget",
                            "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/%s/%s.fasta.gz" % (op.target_database, op.target_database),
                            "-P", fo.database_directory],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["uniprot"]["error"] = True
            sys.stderr.write("Error: There was a problem downloading the UniProt target database in preparation for processing. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
            source_exists = False
    else:
        sys.stdout.write("UniRef90 database exists, skipping download.\n")

    if (not os.path.isfile(op.target_database_file) and source_exists) or (source_exists and args.force):
        sys.stdout.write("Creating UniRef90 DIAMOND database.\n")
        try:
            # Tag output database file with a ".tmp" tag and remove it once the database has been generated
            # just in case something interrupts the database creation process
            subprocess.run([config.steps.map_reads_to_genes.required_programs["diamond"],
                            "makedb",
                            "--in", fo.database_directory + op.target_database + ".fasta.gz",
                            "-d", fo.database_directory + op.target_database + ".tmp"],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
            temp_database_file = fo.database_directory + op.target_database + ".tmp" + op.diamond_db_suffix
            permanent_file = re.sub("\\.tmp", "", temp_database_file)
            subprocess.run(["mv",
                            temp_database_file,
                            permanent_file],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["uniprot"]["error"] = True
            sys.stderr.write("Error: There was a problem creating a DIAMOND database for the  UniProt target database. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
    else:
        sys.stdout.write("UniRef90 DIAMOND database exists, skipping DIAMOND database generation.\n")

    if (not os.path.isfile(op.gene_normalization_file) and source_exists) or (source_exists and args.force):
        sys.stdout.write("Creating UniRef90 gene length table.\n")
        try:
            # Tag output gene length table with a ".tmp" tag and remove it once the table has been generated
            # just in case something interrupts the table creation process
            subprocess.run([fo.source_directory + "create_gene_length_table.py",
                            fo.database_directory + op.target_database + ".fasta.gz",
                            "--output", op.gene_normalization_file + ".tmp"],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
            subprocess.run(["mv",
                            op.gene_normalization_file + ".tmp",
                            op.gene_normalization_file],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["uniprot"]["error"] = True
            sys.stderr.write("Error: There was a problem creating a gene length table for the  UniProt target database. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
    else:
        sys.stdout.write("UniRef90 gene length table exists, skipping gene length table generation.\n")

    source_exists = True
    if not os.path.isfile(fo.database_directory + "idmapping.dat.gz") or args.force:
        sys.stdout.write("Downloading UniRef90 gene-to-ortholog mapping.\n")
        try:
            subprocess.run(["wget",
                            "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz",
                            "-P", fo.database_directory],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["uniprot"]["error"] = True
            sys.stderr.write(
                "Error: There was a problem downloading the UniProt gene-to-ortholog mapping in preparation for processing. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
            source_exists = False
    else:
        sys.stdout.write("UniRef90 gene-to-ortholog mapping exists, skipping download.\n")

    if (not os.path.isfile(op.gene_to_ortholog_file) and source_exists) or (source_exists and args.force):
        sys.stdout.write("Formatting UniRef90 gene-to-ortholog mapping.\n")
        try:
            # Tag output gene length table with a ".tmp" tag and remove it once the table has been generated
            # just in case something interrupts the table creation process
            subprocess.run([fo.source_directory + "create_uniref_gene_to_ortholog.py",
                            fo.database_directory + "idmapping.dat.gz",
                            op.target_database,
                            op.target_ortholog,
                            "--output", op.gene_to_ortholog_file + ".tmp"],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
            subprocess.run(["mv",
                            op.gene_to_ortholog_file + ".tmp",
                            op.gene_to_ortholog_file],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["uniprot"]["error"] = True
            sys.stderr.write("Error: There was a problem formatting the gene-to-ortholog table for the  UniProt target database. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
    else:
        sys.stdout.write("Formatted UniRef90 gene-to-ortholog mapping exists, skipping formatting.\n")
sys.stdout.write("\n")

sys.stdout.write("*" * 80 + "\n")
sys.stdout.write("SETUP COMPLETE\n")
sys.stdout.write("Automated reference data processing has finished. Please see above for details on supporting data files that could not be generated.\n\n")

# Report any database files that are missing
if not os.path.isfile(op.host_index + ".1.bt2"):
    sys.stderr.write("Warning: The host database index (%s) does not exist. You will be unable to perform the default host filtering step without it. Make sure that 'index_directory' in config.file_organization and 'host_database' in config.operation.py are correct.\n" % op.host_index)

if not os.path.isfile(op.target_database_file):
    sys.stderr.write("Warning: The target DIAMOND database (%s) does not exist. You will be unable to perform the default read mapping step without it. Make sure that 'database_directory' in config.file_organization and 'target_database_file' in config.operation.py are correct.\n" % op.target_database_file)

if not os.path.isfile(op.gene_normalization_file):
    sys.stderr.write("Warning: The target database gene normalization table (%s) does not exist. You will be unable to perform the default gene mapping step without it. Make sure that 'gene_normalization_directory' in config.file_organization and 'target_database_file' in config.operation.py are correct.\n" % op.gene_normalization_file)

if not os.path.isfile(op.gene_to_ortholog_file):
    sys.stderr.write("Warning: The gene-to-ortholog table (%s) does not exist. You will be unable to perform the default ortholog mapping step (and the hit filtering step if using either the 'best_ortholog' or 'best_n_orthologs' methods) without it. Make sure that 'gene_to_ortholog_directory' in config.file_organization and 'gene_to_ortholog_file' in config.operation.py are correct.\n" % op.gene_to_ortholog_file)

for ortholog_to_grouping in op.ortholog_to_grouping_files:
    if not os.path.isfile(ortholog_to_grouping):
        sys.stderr.write("Warning: The ortholog-to-grouping map (%s) does not exist. You will be unable to perform the default ortholog aggregation step without it. Make sure that 'ortholog_to_grouping_directory' in config.file_organization and 'ortholog_to_grouping' in config.operation.py are correct.\n" % ortholog_to_grouping)

sys.stdout.write("End of database preparation.\n")
