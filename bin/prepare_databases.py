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
parser.add_argument("--cleanup",
                    "-c",
                    action="store_true",
                    help="Normally, this script will keep all intermediate files to avoid downloading source files and re-processing files in case the running of the script is interrupted or the user decides to change the base configuration that then requires re-processing of the source data (e.g. changing the target UniRef version.). However, the intermediate files can be rather large. Use this option to force the script to remove intermediate files as they are successfully processed, which will save disk space.")

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

sys.stdout.write("Beginning processing of reference data.%s%s" % (os.linesep, os.linesep))
# Try to create and process default database files
if args.human_reference:
    source_exists = True
    if (not os.path.isfile(op.host_database_file) and len(glob.glob(op.host_index + "*.bt2")) == 0) or args.force:
        sys.stdout.write("Downloading human reference.%s" % os.linesep)
        try:
            subprocess.run(["wget",
                            op.host_database_source,
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
            sys.stderr.write("Error: There was a problem downloading the human reference database in preparation for processing. Please see %s and %s for processing logs.%s" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"], os.linesep))
            source_exists = False
    else:
        sys.stdout.write("Human reference exists, skipping download.%s" % os.linesep)

    if source_exists:
        if len(glob.glob(op.host_index + "*.bt2")) == 0 or args.force:
            sys.stdout.write("Creating human reference Bowtie 2 index.%s" % os.linesep)
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
                if args.cleanup:
                    subprocess.run(["rm", op.host_database_file],
                                   stdout=open(processing_results["human_reference"]["stdout"], "a"),
                                   stderr=open(processing_results["human_reference"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["human_reference"]["error"] = True
                sys.stderr.write("Error: There was a problem creating an index for the human reference database. Please see %s and %s for processing logs.%s" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"], os.linesep))
        else:
            sys.stdout.write("Human reference Bowtie 2 index exists, skipping database indexing.%s" % os.linesep)
    sys.stdout.write(os.linesep)

if args.ko_mappings:
    missing_mappings = []
    for ortholog_to_grouping_mapping in op.ortholog_to_grouping_mappings:
        if not os.path.isfile(fo.ortholog_to_grouping_directory +
                              ortholog_to_grouping_mapping +
                              op.ortholog_to_grouping_suffix) or args.force:
            missing_mappings.append(ortholog_to_grouping_mapping)
    if len(missing_mappings) > 0:
        sys.stdout.write("Downloading ortholog-to-grouping mappings.%s" % os.linesep)
        for ortholog_to_grouping_mapping in missing_mappings:
            if not os.path.isfile(fo.ortholog_to_grouping_directory +
                                  ortholog_to_grouping_mapping +
                                  op.ortholog_to_grouping_suffix) or args.force:
                try:
                    subprocess.run(["wget",
                                    op.ortholog_to_grouping_mapping_source + ortholog_to_grouping_mapping,
                                    "-O",
                                    fo.ortholog_to_grouping_directory + ortholog_to_grouping_mapping + op.ortholog_to_grouping_suffix],
                                   stdout=open(processing_results["ortholog_mapping"]["stdout"], "a"),
                                   stderr=open(processing_results["ortholog_mapping"]["stderr"], "a"))
                except (EnvironmentError, subprocess.CalledProcessError):
                    processing_error = True
                    processing_results["ortholog_mapping"]["error"] = True
                    sys.stderr.write("Error: There was a problem downloading the %s ortholog-to-grouping mapping file. Please see %s and %s for processing logs.%s" % (ortholog_to_grouping_mapping, processing_results["ortholog_mapping"]["stdout"], processing_results["ortholog_mapping"]["stderr"], os.linesep))
    else:
        sys.stdout.write("Ortholog-to-grouping mappings exist, skipping download.%s" % os.linesep)
sys.stdout.write(os.linesep)

if args.uniprot:
    source_exists = True
    if (not os.path.isfile(fo.database_directory + "uniref%s.tar.gz" % op.target_database_release) and not os.path.isfile(fo.database_directory + op.target_database + ".tar") and not os.path.isfile(fo.database_directory + op.target_database + ".xml.gz") and not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz") and not os.path.isfile(fo.database_directory + op.target_database + ".uniprot.mapping.gz") and (not os.path.isfile(op.target_database_file) or not os.path.isfile(op.gene_normalization_file) or not os.path.isfile(op.gene_to_ortholog_file))) or args.force:

        sys.stdout.write("Downloading UniRef archive.%s" % os.linesep)
        try:
            subprocess.run(["wget",
                            "%srelease-%s/uniref/uniref%s.tar.gz" % (op.target_database_source, op.target_database_release, op.target_database_release),
                            "-P", fo.database_directory],
                           stdout=open(processing_results["uniprot"]["stdout"], "a"),
                           stderr=open(processing_results["uniprot"]["stderr"], "a"))
        except (EnvironmentError, subprocess.CalledProcessError):
            processing_error = True
            processing_results["uniprot"]["error"] = True
            sys.stderr.write("Error: There was a problem downloading the UniRef archive in preparation for processing. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
            source_exists = False
    else:
        sys.stdout.write("UniRef archive (or later processed files) exists, skipping download.%s" % os.linesep)

    if source_exists:
        if (not os.path.isfile(fo.database_directory + op.target_database + ".tar") and not os.path.isfile(fo.database_directory + op.target_database + ".xml.gz") and not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz") and not os.path.isfile(fo.database_directory + op.target_database + ".uniprot.mapping.gz") and (not os.path.isfile(op.target_database_file) or not os.path.isfile(op.gene_normalization_file) or not os.path.isfile(op.gene_to_ortholog_file))) or args.force:

            sys.stdout.write("Extracting target UniRef database archive from full UniRef archive.%s" % os.linesep)
            try:
                subprocess.run(["tar", "-C", fo.database_directory, "-zxf",
                                fo.database_directory + "uniref%s.tar.gz" % op.target_database_release,
                                op.target_database + ".tar"],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                if args.cleanup:
                    subprocess.run(["rm", fo.database_directory + "uniref%s.tar.gz" % op.target_database_release],
                                   stdout=open(processing_results["uniprot"]["stdout"], "a"),
                                   stderr=open(processing_results["uniprot"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem extracting the target UniRef database archive from the full UniRef archive. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
                source_exists = False
        else:
            sys.stdout.write("Target UniRef database archive (or later processed files) exists, skipping extraction.%s" % os.linesep)

    if source_exists:
        if (not os.path.isfile(fo.database_directory + op.target_database + ".xml.gz") and not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz") and not os.path.isfile(fo.database_directory + op.target_database + ".uniprot.mapping.gz") and (not os.path.isfile(op.target_database_file) or not os.path.isfile(op.gene_normalization_file) or not os.path.isfile(op.gene_to_ortholog_file))) or args.force:

            sys.stdout.write("Extracting target UniRef database XML from target UniRef database archive.%s" % os.linesep)
            try:
                subprocess.run(["tar", "-C", fo.database_directory, "-xf",
                                fo.database_directory + op.target_database + ".tar",
                                op.target_database + ".xml.gz"],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                if args.cleanup:
                    subprocess.run(["rm", fo.database_directory + op.target_database + ".tar"],
                                   stdout=open(processing_results["uniprot"]["stdout"], "a"),
                                   stderr=open(processing_results["uniprot"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem extracting the target UniRef database XML from the target UniRef database archive. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
                source_exists = False
        else:
            sys.stdout.write("Target UniRef database XML (or later processed files) exists, skipping extraction.%s" % os.linesep)

    if source_exists:
        if (not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz") and (not os.path.isfile(op.target_database_file) or not os.path.isfile(op.gene_normalization_file))) or (not os.path.isfile(fo.database_directory + op.target_database + ".uniprot.mapping.gz") and not os.path.isfile(op.gene_to_ortholog_file)) or args.force:

            sys.stdout.write("Parsing target UniRef database XML to create target UniRef FASTA and UniRef-to-UniProt mapping files.%s" % os.linesep)
            try:
                subprocess.run([fo.source_directory + "parse_uniref_xml.py",
                                fo.database_directory + op.target_database + ".xml.gz",
                                fo.database_directory + op.target_database + ".fasta",
                                fo.database_directory + op.target_database + ".uniprot.mapping"],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                subprocess.run(["gzip",
                                fo.database_directory + op.target_database + ".fasta",
                                fo.database_directory + op.target_database + ".uniprot.mapping"],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                if args.cleanup:
                    subprocess.run(["rm", fo.database_directory + op.target_database + ".xml.gz"],
                                   stdout=open(processing_results["uniprot"]["stdout"], "a"),
                                   stderr=open(processing_results["uniprot"]["stderr"], "a"))
                    if os.path.isfile(op.target_database_file) and os.path.isfile(op.gene_normalization_file):
                        subprocess.run(["rm", fo.database_directory + op.target_database + ".fasta.gz"])
                    if os.path.isfile(op.gene_to_ortholog_file):
                        subprocess.run(["rm", fo.database_directory + op.target_database + ".uniprot.mapping.gz"])
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem parsing the target UniRef database XML. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
                source_exists = False
        else:
            sys.stdout.write("Target UniRef FASTA and UniRef-to-UniProt mapping files (or later processed files) exist, skipping parsing.%s" % os.linesep)

    if source_exists:
        if not os.path.isfile(op.target_database_file) or args.force:
            sys.stdout.write("Creating target UniRef DIAMOND database.%s" % os.linesep)
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
                sys.stderr.write("Error: There was a problem creating a DIAMOND database for the target UniProt database. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
        else:
            sys.stdout.write("Target UniRef DIAMOND database exists, skipping DIAMOND database generation.%s" % os.linesep)

    if source_exists:
        if not os.path.isfile(op.gene_normalization_file) or args.force:
            sys.stdout.write("Creating target UniRef gene length table.%s" % os.linesep)
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
                if args.cleanup and os.path.isfile(op.target_database_file):
                    subprocess.run(["rm", fo.database_directory + op.target_database + ".fasta.gz"],
                                   stdout=open(processing_results["uniprot"]["stdout"], "a"),
                                   stderr=open(processing_results["uniprot"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem creating a gene length table for the target UniProt database. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
        else:
            sys.stdout.write("Target UniRef gene length table exists, skipping gene length table generation.%s" % os.linesep)

    if source_exists:
        if (not os.path.isfile(fo.database_directory + "uniprot_sprot-only%s.tar.gz" % op.target_database_release) and not os.path.isfile(fo.database_directory + "uniprot_sprot.dat.gz") and not os.path.isfile(op.gene_to_ortholog_file)) or args.force:
            sys.stdout.write("Downloading UniProt mapping archive.%s" % os.linesep)
            try:
                subprocess.run(["wget",
                                "%s/release-%s/knowledgebase/uniprot_sprot-only%s.tar.gz" % (op.target_database_source, op.target_database_release, op.target_database_release),
                                "-P", fo.database_directory],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write(
                    "Error: There was a problem downloading the UniProt mapping archive in preparation for processing. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
                source_exists = False
        else:
            sys.stdout.write("UniProt mapping archive (or later processed files) exists, skipping download.%s" % os.linesep)

    if source_exists:
        if (not os.path.isfile(fo.database_directory + "uniprot_sprot.dat.gz") and not os.path.isfile(op.gene_to_ortholog_file)) or args.force:
            sys.stdout.write("Extracting UniProt gene-to-ortholog mapping.%s" % os.linesep)
            try:
                subprocess.run(["tar",
                                "-C", fo.database_directory,
                                "-zxf", fo.database_directory + "uniprot_sprot-only%s.tar.gz" % op.target_database_release,
                                "uniprot_sprot.dat.gz"],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                if args.cleanup:
                    subprocess.run(["rm", fo.database_directory + "uniprot_sprot-only%s.tar.gz" % op.target_database_release],
                                   stdout=open(processing_results["uniprot"]["stdout"], "a"),
                                   stderr=open(processing_results["uniprot"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write(
                    "Error: There was a problem extracting the UniProt gene-to-ortholog mapping in preparation for processing. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
                source_exists = False
        else:
            sys.stdout.write("UniProt gene-to-ortholog mapping (or later processed files) exists, skipping download.%s" % os.linesep)

    if source_exists:
        if not os.path.isfile(op.gene_to_ortholog_file) or args.force:
            sys.stdout.write("Formatting target UniRef gene-to-ortholog mapping.%s" % os.linesep)
            try:
                # Tag output gene length table with a ".tmp" tag and remove it once the table has been generated
                # just in case something interrupts the table creation process
                subprocess.run([fo.source_directory + "create_uniref_gene_to_ortholog.py",
                                fo.database_directory + op.target_database + ".uniprot.mapping.gz",
                                fo.database_directory + "uniprot_sprot.dat.gz",
                                op.target_ortholog,
                                "--output", op.gene_to_ortholog_file + ".tmp"],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                subprocess.run(["mv",
                                op.gene_to_ortholog_file + ".tmp",
                                op.gene_to_ortholog_file],
                               stdout=open(processing_results["uniprot"]["stdout"], "a"),
                               stderr=open(processing_results["uniprot"]["stderr"], "a"))
                if args.cleanup:
                    subprocess.run(["rm",
                                    fo.database_directory + op.target_database + ".uniprot.mapping.gz",
                                    fo.database_directory + "uniprot_sprot.dat.gz"],
                                   stdout=open(processing_results["uniprot"]["stdout"], "a"),
                                   stderr=open(processing_results["uniprot"]["stderr"], "a"))
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem formatting the gene-to-ortholog table for the target UniRef database. Please see %s and %s for processing logs.%s" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"], os.linesep))
        else:
            sys.stdout.write("Formatted target UniRef gene-to-ortholog mapping exists, skipping formatting.%s" % os.linesep)
    sys.stdout.write(os.linesep)

sys.stdout.write("*" * 80 + os.linesep)
sys.stdout.write("SETUP COMPLETE%s" % os.linesep)
sys.stdout.write("Automated reference data processing has finished. Please see above for details on supporting data files that could not be generated.%s%s" % (os.linesep, os.linesep))

# Report any database files that are missing
if not os.path.isfile(op.host_index + ".1.bt2"):
    sys.stderr.write("Warning: The host database index (%s) does not exist. You will be unable to perform the default host filtering step without it. Make sure that 'index_directory' in config.file_organization and 'host_database' in config.operation.py are correct.%s" % (op.host_index, os.linesep))

if not os.path.isfile(op.target_database_file):
    sys.stderr.write("Warning: The target DIAMOND database (%s) does not exist. You will be unable to perform the default read mapping step without it. Make sure that 'database_directory' in config.file_organization and 'target_database_file' in config.operation.py are correct.%s" % (op.target_database_file, os.linesep))

if not os.path.isfile(op.gene_normalization_file):
    sys.stderr.write("Warning: The target database gene normalization table (%s) does not exist. You will be unable to perform the default gene mapping step without it. Make sure that 'gene_normalization_directory' in config.file_organization and 'target_database_file' in config.operation.py are correct.%s" % (op.gene_normalization_file, os.linesep))

if not os.path.isfile(op.gene_to_ortholog_file):
    sys.stderr.write("Warning: The gene-to-ortholog table (%s) does not exist. You will be unable to perform the default ortholog mapping step (and the hit filtering step if using either the 'best_ortholog' or 'best_n_orthologs' methods) without it. Make sure that 'gene_to_ortholog_directory' in config.file_organization and 'gene_to_ortholog_file' in config.operation.py are correct.%s" % (op.gene_to_ortholog_file, os.linesep))

for ortholog_to_grouping in op.ortholog_to_grouping_files:
    if not os.path.isfile(ortholog_to_grouping):
        sys.stderr.write("Warning: The ortholog-to-grouping map (%s) does not exist. You will be unable to perform the default ortholog aggregation step without it. Make sure that 'ortholog_to_grouping_directory' in config.file_organization and 'ortholog_to_grouping' in config.operation.py are correct.%s" % (ortholog_to_grouping, os.linesep))

sys.stdout.write("End of database preparation.%s" % os.linesep)
