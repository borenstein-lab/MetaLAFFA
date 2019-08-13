import argparse
import sys
import os
import subprocess
from config import env
import config.file_organization as fo
import config.operation as op
import config.steps
import pkgutil
import importlib

parser = argparse.ArgumentParser(description="Setup script for installing third party tools and downloading default databases.")
parser.add_argument("--pip", "-p", default="pip3", help="Which pip to run for local python package installation (default: %(default)s)")
parser.add_argument("--no_snakemake", "-ns", action="store_true", help="If used, do not install Snakemake locally. Note that Snakemake is required for pipeline operation, so only use this option if you have already installed Snakemake.")
parser.add_argument("--no_blast", "-nbl", action="store_true", help="If used, do not download and install it in the pipeline directory for host filtering.")
parser.add_argument("--no_bmtagger", "-nbm", action="store_true", help="If used, do not download BMTagger and install it in the pipeline directory for host filtering.")
parser.add_argument("--no_human_reference", "-nh", action="store_true", help="If used, do not download the human reference provided with BMTagger (hs37).")
parser.add_argument("--no_markduplicates", "-nma", action="store_true", help="If used, do not download PICARD and samtools and install them in the pipeline directory for duplicate filtering.")
parser.add_argument("--no_trimmomatic", "-nt", action="store_true", help="If used, do not download Trimmomatic and install it in the pipeline directory for quality filtering.")
parser.add_argument("--no_diamond", "-nd", action="store_true", help="If used, do not download DIAMOND and install it in the pipeline directory for read mapping.")
parser.add_argument("--no_musicc", "-nmu", action="store_true", help="If used, do not install MUSiCC via pip for ortholog abundance correction.")
parser.add_argument("--no_ko_mappings", "-nkm", action="store_true", help="If used, to not download bacterial ko-to-module and ko-to-pathway mappings from the 2013 version of the KEGG database.")
parser.add_argument("--uniprot", "-u", action="store_true", help="If used, download the default UNIPROT gene sequence reference database for read mapping.")
parser.add_argument("--no_jobscript", "-nj", action="store_true", help="If used, do not configure the template jobscript to standardize the environment when running cluster jobs remotely")

args = parser.parse_args()

# Create any missing directories
for required_directory in fo.required_directories:
    if not os.path.isdir(required_directory):
        os.makedirs(required_directory)

install_error = False
pip_installed = True
make_installed = True
cmake_installed = True
installed_tool_dependencies = {
    "bmtools": ["human reference bitmask"],
    "srprism": ["human reference index"],
    "blast": ["human reference blast database"],
    "diamond": ["annotated gene diamond database"]
}
install_results = {}
for tool in ["python_package", "blast", "bmtools", "srprism", "samtools", "diamond"]:
    install_results[tool] = {}
    install_results[tool]["error"] = False
    install_results[tool]["stdout"] = fo.source_directory + tool + "_install.out"
    install_results[tool]["stderr"] = fo.source_directory + tool + "_install.err"

if subprocess.run(["which", args.pip], capture_output=True, env=env).returncode != 0:
    pip_installed = False

if subprocess.run(["which", "make"], capture_output=True, env=env).returncode != 0:
    make_installed = False

if subprocess.run(["which", "cmake"], capture_output=True, env=env).returncode != 0:
    cmake_installed = False

# Report any missing required tools
if not pip_installed:
    sys.stderr.write(args.pip + " not found. pip is required for installation of [snakemake, musicc]. These tools will not be installed during setup.\n")

if not make_installed:
    sys.stderr.write("make not found. Make is required for installation of [blast, bmtagger, markduplicates, diamond]. These tools will not be installed during setup.\n")

if not cmake_installed:
    sys.stderr.write("cmake not found. CMake is required for installation of [diamond]. These tools will not be installed during setup.\n")

# Try to install third party tools
if pip_installed:
    pip_install_list = []
    if not args.no_snakemake:
        pip_install_list.append("snakemake")
        pip_install_list.append("psutil")
    if not args.no_musicc:
        pip_install_list.append("musicc")
    try:
        subprocess.run([args.pip, "install", "-t", fo.python_package_directory] + pip_install_list, stdout=open(install_results["python_package"]["stdout"], "w"), stderr=open(install_results["python_package"]["stderr"], "w"), env=env)
    except (EnvironmentError, subprocess.CalledProcessError):
        install_error = True
        install_results["python_package"]["error"] = True

if not args.no_blast and make_installed:
    if not os.path.isdir(fo.source_directory + "ncbi-blast-2.2.31+-src"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-src.tar.gz", "-P", fo.source_directory], env=env)
        subprocess.run(["tar", "-zxf", "ncbi-blast-2.2.31+-src.tar.gz"], cwd=fo.source_directory, env=env)

    if not os.path.isdir(config.steps.host_filter.required_programs["blastn_dir"]) or not os.path.isfile(config.steps.host_filter.required_programs["blastn_dir"] + "makeblastdb"):
        subprocess.run(["./configure"], cwd=fo.source_directory + "ncbi-blast-2.2.31+-src/c++/", env=env)

        try:
            subprocess.run(["make", "all_r"], cwd=fo.source_directory + "ncbi-blast-2.2.31+-src/c++/ReleaseMT/build/", stdout=open(install_results["blast"]["stdout"], "w"), stderr=open(install_results["blast"]["stderr"], "w"), env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["blast"]["error"] = True

if not args.no_bmtagger and make_installed:
    if not os.path.isdir(fo.source_directory + "bmtools/"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmtools.tar.gz", "-P", fo.source_directory], env=env)
        subprocess.run(["tar", "-zxf", "bmtools.tar.gz"], cwd=fo.source_directory, env=env)
    if not os.path.isfile(config.steps.host_filter.required_programs["bmtool"]) or not os.path.isfile(config.steps.host_filter.required_programs["bmfilter"]) or not os.path.isfile(config.steps.host_filter.required_programs["bmtagger"]) or not os.path.isfile(config.steps.host_filter.required_programs["extract_fa"]):

        try:
            subprocess.run(["make"], cwd=fo.source_directory + "bmtools/", stdout=open(install_results["bmtools"]["stdout"], "w"), stderr=open(install_results["bmtools"]["stderr"], "w"), env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["bmtools"]["error"] = True

    if not os.path.isdir(fo.source_directory + "srprism/"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/src/srprism.tar.gz", "-P", fo.source_directory], env=env)
        subprocess.run(["tar", "-zxf", "srprism.tar.gz"], cwd=fo.source_directory, env=env)
    if not os.path.isfile(config.steps.host_filter.required_programs["srprism"]):
        subprocess.run(["./ac.sh"], cwd=fo.source_directory + "srprism/gnuac/", env=env)
        subprocess.run(["./configure"], cwd=fo.source_directory + "srprism/gnuac/", env=env)

        try:
            subprocess.run(["make"], cwd=fo.source_directory + "srprism/gnuac/", stdout=open(install_results["srprism"]["stdout"], "w"), stderr=open(install_results["srprism"]["stderr"], "w"), env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["srprism"]["error"] = True

if not args.no_markduplicates:
    if not os.path.isfile(fo.source_directory + "picard.jar"):
        subprocess.run(["wget", "https://github.com/broadinstitute/picard/releases/download/2.20.2/picard.jar", "-P", fo.source_directory], env=env)

    if not os.path.isdir(fo.source_directory + "samtools-1.9/"):
        subprocess.run(["wget", "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2", "-P", fo.source_directory], env=env)
        subprocess.run(["tar", "-jxf", "samtools-1.9.tar.bz2"], cwd=fo.source_directory, env=env)

    if not os.path.isfile(config.steps.duplicate_filter.required_programs["samtools"]) and make_installed:
        subprocess.run(["./configure", "--without-curses", "--disable-lzma"], cwd=fo.source_directory + "samtools-1.9/", env=env)
        try:
            subprocess.run(["make"], cwd=fo.source_directory + "samtools-1.9/", stdout=open(install_results["samtools"]["stdout"], "w"), stderr=open(install_results["samtools"]["stderr"], "w"), env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["samtools"]["error"] = True

if not args.no_trimmomatic:
    if not os.path.isfile(config.steps.quality_filter.required_programs["trimmer"]):
        subprocess.run(["wget", "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip", "-P", fo.source_directory], env=env)
        subprocess.run(["unzip", fo.source_directory + "Trimmomatic-0.39.zip", "-d", fo.source_directory], env=env)

if not args.no_diamond and make_installed and cmake_installed:
    if not os.path.isdir(fo.source_directory + "diamond-0.9.22/"):
        subprocess.run(["wget", "http://github.com/bbuchfink/diamond/archive/v0.9.22.tar.gz", "-P", fo.source_directory], env=env)
        subprocess.run(["tar", "-zxf", "v0.9.22.tar.gz"], cwd=fo.source_directory, env=env)

    if not os.path.isfile(config.steps.map_reads_to_genes.required_programs["diamond"]):
        try:
            subprocess.run(["cmake", ".", "-DCMAKE_INSTALL_PREFIX=" + os.getcwd() + "/" + fo.source_directory + "diamond-0.9.22/"], cwd=fo.source_directory + "diamond-0.9.22/", stdout=open(install_results["diamond"]["stdout"], "w"), stderr=open(install_results["diamond"]["stderr"], "w"), env=env)
            subprocess.run(["make", "install"], cwd=fo.source_directory + "diamond-0.9.22/", stdout=open(install_results["diamond"]["stdout"], "a"), stderr=open(install_results["diamond"]["stderr"], "a"), env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["diamond"]["error"] = True

# Report any installation errors
if install_error:
    for key in install_results.keys():
        if install_results[key]["error"]:
            sys.stderr.write("%s failed to install properly. See %s and %s for installation details.\n" % (key, install_results[key]["stdout"], install_results[key]["stderr"]))
            if len(installed_tool_dependencies[key]) > 0:
                sys.stderr.write("The following were not generated because %s is required for their generation: %s\n" % (key, ", ".join(installed_tool_dependencies[key])))
            sys.stderr.write("\n")

# Report any tools that are missing
all_missing_programs = set()
if subprocess.run(["which", op.python], capture_output=True, env=env).returncode != 0:
    sys.stderr.write("Warning: %s is not present in your current environment. This will be required for running Python scripts. Please either install Python or change the 'python' variable in the config.operation.py submodule to point to a valid Python executable.\n" % op.python)
    all_missing_programs.add(op.python)

if subprocess.run(["which", op.java], capture_output=True, env=env).returncode != 0:
    sys.stderr.write("Warning: %s is not present in your current environment. This will be required for running Java jars. Please either install Java or change the 'java' variable in the config.operation.py submodule to point to a valid Java executable.\n" % op.java)
    all_missing_programs.add(op.java)

if subprocess.run(["which", op.snakemake], capture_output=True, env=env).returncode != 0:
    sys.stderr.write("Warning: %s is not present in your current environment. This will be required for running metaLAFFA. Please either install Snakemake or change the 'snakemake' variable in the config.operation.py submodule to point to a valid snakemake executable.\n" % op.snakemake)
    all_missing_programs.add(op.snakemake)

makeblastdb_location = config.steps.host_filter.required_programs["blastn_dir"] + "makeblastdb"
if subprocess.run(["which", makeblastdb_location], capture_output=True, env=env).returncode != 0:
    sys.stderr.write("Warning: %s is not present and executable in your current environment. This will be required for processing the default human reference database for host filtering. Check for BLAST+ installation errors if you tried to install BLAST+ with the setup script.\n" % makeblastdb_location)
    all_missing_programs.add(makeblastdb_location)

for importer, modname, ispkg in pkgutil.iter_modules(config.steps.__path__):
    required_programs = importlib.import_module(".".join(["config.steps", modname])).required_programs.values()

    # Check if each file exists or is an executable in the environment
    step_missing_programs = []
    for required_program in required_programs:
        program_exists = os.path.exists(required_program)
        executable_exists = subprocess.run(["which", required_program], capture_output=True, env=env).returncode == 0
        if not program_exists and not executable_exists:
            step_missing_programs.append(required_program)

    # If a program is missing, inform the user and add it to the set of missing programs
    if len(step_missing_programs) > 0:
        sys.stderr.write("Warning: program(s) missing for the %s pipeline step. The step may not run correctly without the program(s):\n" % modname)
        for missing_program in step_missing_programs:
            sys.stderr.write("%s\n" % missing_program)
        all_missing_programs.union(set(step_missing_programs))

# Try to create and process default database files
if not args.no_human_reference and config.steps.host_filter.required_programs["bmtool"] not in all_missing_programs and config.steps.host_filter.required_programs["srprism"] not in all_missing_programs and makeblastdb_location not in all_missing_programs:
    if not os.path.isfile(fo.database_directory + op.host_database + ".fa"):
        subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/%s.fa" % op.host_database, "-P", fo.database_directory], env=env)

    if not os.path.isfile(op.host_bitmask_file) and not install_results["bmtools"]["error"]:
        subprocess.run([config.steps.host_filter.required_programs["bmtool"], "-d", op.host_database_file + ".fa", "-o", op.host_bitmask_file, "-A", "0", "-w", "18"], env=env)

    if not os.path.isfile(op.host_index_file + ".idx") and not install_results["srprism"]["error"]:
        subprocess.run([config.steps.host_filter.required_programs["srprism"], "mkindex", "-i", op.host_database_file + ".fa", "-o", op.host_index_file, "-M", "7168"], env=env)

    if not os.path.isfile(op.host_database_file + op.blast_db_suffix) and not install_results["blast"]["error"]:
        subprocess.run([makeblastdb_location, "-in", op.host_database_file + ".fa", "-dbtype", "nucl", "-out", op.host_database_file], env=env)

if not args.no_ko_mappings:
    for ortholog_to_grouping_mapping in op.ortholog_to_grouping_mappings:
        if not os.path.isfile(fo.ortholog_to_grouping_directory + ortholog_to_grouping_mapping + op.ortholog_to_grouping_suffix):
            subprocess.run(["wget", "https://github.com/borenstein-lab/fishtaco/raw/master/fishtaco/data/" + ortholog_to_grouping_mapping, "-O", fo.ortholog_to_grouping_directory + ortholog_to_grouping_mapping + op.ortholog_to_grouping_suffix], env=env)

if args.uniprot and config.steps.map_reads_to_genes.required_programs["diamond"] not in all_missing_programs and op.python not in all_missing_programs:
    if not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz"):
        subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/%s/%s.fasta.gz" % (op.target_database, op.target_database), "-P", fo.database_directory], env=env)

    if not os.path.isfile(op.target_database_file) and not install_results["diamond"]["error"]:
        subprocess.run([config.steps.map_reads_to_genes.required_programs["diamond"], "makedb", "--in", fo.database_directory + op.target_database + ".fasta.gz", "-d", fo.database_directory + op.target_database], env=env)

    if not os.path.isfile(fo.database_directory + "idmapping.dat.gz"):
        subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz", "-P", fo.database_directory], env=env)

    if not os.path.isfile(op.gene_normalization_file):
        subprocess.run([op.python, fo.source_directory + "create_gene_length_table.py", fo.database_directory + op.target_database + ".fasta.gz", "--output", op.gene_normalization_file], env=env)

    if not os.path.isfile(op.gene_to_ortholog_file):
        subprocess.run([op.python, fo.source_directory + "create_uniref_gene_to_ortholog.py", fo.database_directory + "idmapping.dat.gz", op.target_database, op.target_ortholog, "--output", op.gene_to_ortholog_file], env=env)

# Report any database files that are missing
if not os.path.isfile(op.host_bitmask_file):
    sys.stderr.write("Warning: The host database bitmask (%s) does not exist. You will be unable to perform the default host filtering step without it. Make sure that 'bitmask_directory' in config.file_organization and 'host_database' in config.operation.py are correct.\n" % op.host_bitmask_file)

if not os.path.isfile(op.host_database_file + op.blast_db_suffix):
    sys.stderr.write("Warning: The host BLAST database sequence file (%s) does not exist. You will be unable to perform the default host filtering step without it. Make sure that 'database_directory' in config.file_organization and 'host_database' in config.operation.py are correct.\n" % op.host_database_file + op.blast_db_suffix)

if not os.path.isfile(op.host_index_file + op.srprism_index_suffix):
    sys.stderr.write("Warning: The host database index (%s) does not exist. You will be unable to perform the default host filtering step without it. Make sure that 'index_directory' in config.file_organization and 'host_database' in config.operation.py are correct.\n" % op.host_index_file)

if not os.path.isfile(op.target_database_file):
    sys.stderr.write("Warning: The target DIAMOND database (%s) does not exist. You will be unable to perform the default read mapping step without it. Make sure that 'database_directory' in config.file_organization and 'target_database_file' in config.operation.py are correct.\n" % op.target_database_file)

if not os.path.isfile(op.gene_normalization_file):
    sys.stderr.write("Warning: The target database gene normalization table (%s) does not exist. You will be unable to perform the default gene mapping step without it. Make sure that 'gene_normalization_directory' in config.file_organization and 'target_database_file' in config.operation.py are correct.\n" % op.gene_normalization_file)

if not os.path.isfile(op.gene_to_ortholog_file):
    sys.stderr.write("Warning: The gene-to-ortholog table (%s) does not exist. You will be unable to perform the default ortholog mapping step (and the hit filtering step if using either the 'best_ortholog' or 'best_n_orthologs' methods) without it. Make sure that 'gene_to_ortholog_directory' in config.file_organization and 'gene_to_ortholog_file' in config.operation.py are correct.\n" % op.gene_to_ortholog_file)

for ortholog_to_grouping in op.ortholog_to_grouping_files:
    if not os.path.isfile(ortholog_to_grouping):
        sys.stderr.write("Warning: The ortholog-to-grouping map (%s) does not exist. You will be unable to perform the default ortholog aggregation step without it. Make sure that 'ortholog_to_grouping_directory' in config.file_organization and 'ortholog_to_grouping' in config.operation.py are correct.\n" % ortholog_to_grouping)

# Configure the job script for running Snakemake on a cluster
if not args.no_jobscript:
    with open(fo.source_directory + "template_jobscript.sh") as template, open(fo.source_directory + "configured_jobscript.sh", "w") as configured:
        for line in template:
            configured.write(line.format(python=op.python, metaLAFFA_directory=os.getcwd(), PYTHONPATH="{{PYTHONPATH}}"))

