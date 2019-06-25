import argparse
import sys
import os
import subprocess
import config.file_organization as fo
import config.operation as op
import config.steps.host_filter
import config.steps.duplicate_filter
import config.steps.quality_filter
import config.steps.map_reads_to_genes

parser = argparse.ArgumentParser(description="Setup script for installing third party tools and downloading default databases.")
parser.add_argument("--no_snakemake", "-ns", action="store_true", help="If used, do not install Snakemake locally. Note that Snakemake is required for pipeline operation, so only use this option if you have already installed Snakemake.")
parser.add_argument("--no_blast", "-nbl", action="store_true", help="If used, do not download and install it in the pipeline directory for host filtering.")
parser.add_argument("--no_bmtagger", "-nbm", action="store_true", help="If used, do not download BMTagger and install it in the pipeline directory for host filtering.")
parser.add_argument("--no_human_reference", "-nh", action="store_true", help="If used, do not download the human reference provided with BMTagger (hs37).")
parser.add_argument("--no_markduplicates", "-nma", action="store_true", help="If used, do not download PICARD and samtools and install them in the pipeline directory for duplicate filtering.")
parser.add_argument("--no_trimmomatic", "-nt", action="store_true", help="If used, do not download Trimmomatic and install it in the pipeline directory for quality filtering.")
parser.add_argument("--no_diamond", "-nd", action="store_true", help="If used, do not download DIAMOND and install it in the pipeline directory for read mapping.")
parser.add_argument("--no_uniprot", "-nu", action="store_true", help="If used, do not download the default UNIPROT gene sequence reference database for read mapping.")
parser.add_argument("--no_musicc", "-nmu", action="store_true", help="If used, do not install MUSiCC via pip for ortholog abundance correction.")

args = parser.parse_args()

# Create any directories required for setup steps
if not args.no_blast or not args.no_bmtagger or not args.no_markduplicates or not args.no_trimmomatic or not args.no_diamond:
    if not os.path.isdir(fo.source_directory):
        os.makedirs(fo.source_directory)

if not args.no_human_reference or not args.no_uniprot:
    if not os.path.isdir(fo.database_directory):
        os.makedirs(fo.database_directory)

if not args.no_human_reference:
    if not os.path.isdir(fo.bitmask_directory):
        os.makedirs(fo.bitmask_directory)
    if not os.path.isdir(fo.srprism_directory):
        os.makedirs(fo.srprism_directory)

if not args.no_snakemake:
    subprocess.run(["pip3", "install", "--user", "snakemake"])

if not args.no_blast:
    if not os.path.isdir(fo.source_directory + "ncbi-blast-2.2.31+-src"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-src.tar.gz", "-P", fo.source_directory])
        subprocess.run(["tar", "-zxf", "ncbi-blast-2.2.31+-src.tar.gz"], cwd=fo.source_directory)

    if not os.path.isdir(config.steps.host_filter.resource_params["blastn_dir"]):
        subprocess.run(["./configure"], cwd=fo.source_directory + "ncbi-blast-2.2.31+-src/c++/")
        subprocess.run(["make", "all_r"], cwd=fo.source_directory + "ncbi-blast-2.2.31+-src/c++/ReleaseMT/build/")


if not args.no_bmtagger:
    if not os.path.isdir(fo.source_directory + "bmtools/"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmtools.tar.gz", "-P", fo.source_directory])
        subprocess.run(["tar", "-zxf", "bmtools.tar.gz"], cwd=fo.source_directory)
    if not os.path.isfile(config.steps.host_filter.resource_params["bmtool"]) or not os.path.isfile(config.steps.host_filter.resource_params["bmfilter"]) or not os.path.isfile(config.steps.host_filter.resource_params["bmtagger"]) or not os.path.isfile(config.steps.host_filter.resource_params["extract_fa"]):
        subprocess.run(["make"], cwd=fo.source_directory + "bmtools/")

    if not os.path.isdir(fo.source_directory + "srpsim/"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/src/srprism.tar.gz", "-P", fo.source_directory])
        subprocess.run(["tar", "-zxf", "srprism.tar.gz"], cwd=fo.source_directory)
    if not os.path.isfile(config.steps.host_filter.resource_params["srprism"]):
        subprocess.run(["./ac.sh"], cwd=fo.source_directory + "srprism/gnuac/")
        subprocess.run(["./configure"], cwd=fo.source_directory + "srprism/gnuac/")
        subprocess.run(["make"], cwd=fo.source_directory + "srprism/gnuac/")

if not args.no_human_reference:
    if not os.path.isfile(fo.database_directory + "hs37.fa"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/hs37.fa", "-P", fo.database_directory])

    if not os.path.isfile(fo.bitmask_directory + "hs37.bitmask"):
        subprocess.run([config.steps.host_filter.resource_params["bmtool"], "-d", fo.database_directory + "hs37.fa", "-o", fo.bitmask_directory + "hs37.bitmask", "-A", "0", "-w", "18"])

    if not os.path.isfile(fo.srprism_directory + "hs37.srprism.idx"):
        subprocess.run([config.steps.host_filter.resource_params["srprism"], "mkindex", "-i", fo.database_directory + "hs37.fa", "-o", fo.srprism_directory + "hs37.srprism", "-M", "7168"])

    if not os.path.isfile(fo.database_directory + "hs37.nhr"):
        subprocess.run([config.steps.host_filter.resource_params["blastn_dir"] + "makeblastdb", "-in", fo.database_directory + "hs37.fa", "-dbtype", "nucl", "-out", fo.database_directory + "hs37"])

if not args.no_markduplicates:
    if not os.path.isfile(fo.source_directory + "picard.jar"):
        subprocess.run(["wget", "https://github.com/broadinstitute/picard/releases/download/2.20.2/picard.jar", "-P", fo.source_directory])

    if not os.path.isdir(fo.source_directory + "samtools-1.9/"):
        subprocess.run(["wget", "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2", "-P", fo.source_directory])
        subprocess.run(["tar", "-jxf", "samtools-1.9.tar.bz2"], cwd=fo.source_directory)

    if not os.path.isfile(config.steps.duplicate_filter.resource_params["samtools"]):
        subprocess.run(["./configure", "--without-curses", "--disable-lzma"], cwd=fo.source_directory + "samtools-1.9/")
        subprocess.run(["make"], cwd=fo.source_directory + "samtools-1.9/")

if not args.no_trimmomatic:
    if not os.path.isfile(config.steps.quality_filter.resource_params["trimmer"]):
        subprocess.run(["wget", "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip", "-P", fo.source_directory])
        subprocess.run(["unzip", fo.source_directory + "Trimmomatic-0.39.zip", "-d", fo.source_directory])

if not args.no_diamond:
    if not os.path.isdir(fo.source_directory + "diamond-0.9.22/"):
        subprocess.run(["wget", "http://github.com/bbuchfink/diamond/archive/v0.9.22.tar.gz", "-P", fo.source_directory])
        subprocess.run(["tar", "-zxf", "v0.9.22.tar.gz"], cwd=fo.source_directory)

    if not os.path.isfile(config.steps.map_reads_to_genes.resource_params["diamond"]):
        subprocess.run(["cmake", ".", "-DCMAKE_INSTALL_PREFIX=" + os.getcwd() + fo.source_directory + "diamond-0.9.22/"], cwd=fo.source_directory + "diamond-0.9.22/")
        subprocess.run(["make", "install"], cwd=fo.source_directory + "diamond-0.9.22/")

if not args.no_uniprot:
    if not os.path.isfile(fo.database_directory + "uniref90.fasta.gz"):
        subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz", "-P", fo.database_directory])

    if not os.path.isfile(fo.database_directory + "uniref90.dmnd"):
        subprocess.run([config.steps.map_reads_to_genes.resource_params["diamond"], "makedb", "--in", fo.database_directory + "uniref90.fasta.gz", "-d", fo.database_directory + "uniref90"])

    if not os.path.isfile(fo.database_directory + "idmapping.dat.gz"):
        subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz", "-P", fo.database_directory])

    if not os.path.isfile(fo.gene_normalization_directory + "uniref90.gene_normalization"):
        subprocess.run([fo.source_directory + "create_gene_length_table.py", fo.database_directory + "uniref90.fasta.gz", "--output", fo.database_directory + "uniref90.gene_normalization"])

    if not os.path.isfile(fo.gene_to_ortholog_directory + "uniref_to_ko"):
        subprocess.run([fo.source_directory + "create_uniref_gene_to_ortholog.py", fo.database_directory + "idmapping.dat.gz", "KO", "--output", fo.gene_to_ortholog_directory + "uniref_to_ko"])

if not args.no_musicc:
    subprocess.run(["pip", "install", "--user", "numpy", "scipy", "scikit-learn==0.17.1", "pandas", "musicc"])
