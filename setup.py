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
parser.add_argument("--pip",
                    "-p",
                    default="pip3",
                    help="Which pip to run for local python package installation (default: %(default)s)")
parser.add_argument("--no_snakemake",
                    "-ns",
                    action="store_true",
                    help="If used, do not install Snakemake locally. Note that Snakemake is required for pipeline operation, so only use this option if you have already installed Snakemake.")
parser.add_argument("--no_blast",
                    "-nbl",
                    action="store_true",
                    help="If used, do not download and install it in the pipeline directory for host filtering.")
parser.add_argument("--no_bmtagger",
                    "-nbm",
                    action="store_true",
                    help="If used, do not download BMTagger and install it in the pipeline directory for host filtering.")
parser.add_argument("--no_human_reference",
                    "-nh",
                    action="store_true",
                    help="If used, do not download the human reference provided with BMTagger (hs37).")
parser.add_argument("--no_markduplicates",
                    "-nma",
                    action="store_true",
                    help="If used, do not download PICARD and samtools and install them in the pipeline directory for duplicate filtering.")
parser.add_argument("--no_trimmomatic",
                    "-nt",
                    action="store_true",
                    help="If used, do not download Trimmomatic and install it in the pipeline directory for quality filtering.")
parser.add_argument("--no_diamond",
                    "-nd",
                    action="store_true",
                    help="If used, do not download DIAMOND and install it in the pipeline directory for read mapping.")
parser.add_argument("--no_musicc",
                    "-nmu",
                    action="store_true",
                    help="If used, do not install MUSiCC via pip for ortholog abundance correction.")
parser.add_argument("--no_empanada",
                    "-ne",
                    action="store_true",
                    help="If used, do not install EMPANADA via pip for ortholog aggregation.")
parser.add_argument("--no_ko_mappings",
                    "-nkm",
                    action="store_true",
                    help="If used, to not download bacterial ko-to-module and ko-to-pathway mappings from the 2013 version of the KEGG database.")
parser.add_argument("--uniprot",
                    "-u",
                    action="store_true",
                    help="If used, download the default UNIPROT gene sequence reference database for read mapping.")
parser.add_argument("--no_jobscript",
                    "-nj",
                    action="store_true",
                    help="If used, do not configure the template jobscript to standardize the environment when running cluster jobs remotely")

args = parser.parse_args()

# Create any missing directories
for required_directory in fo.required_directories:
    if not os.path.isdir(required_directory):
        os.makedirs(required_directory)

install_error = False
processing_error = False
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
for tool in ["python_package", "blast", "bmtools", "srprism", "picard", "samtools", "trimmomatic", "diamond"]:
    install_results[tool] = {}
    install_results[tool]["error"] = False
    install_results[tool]["stdout"] = fo.source_directory + tool + "_install.out"
    install_results[tool]["stderr"] = fo.source_directory + tool + "_install.err"

processing_results = {}
for database in ["human_reference", "ortholog_mapping", "uniprot"]:
    processing_results[database] = {}
    processing_results[database]["error"] = False
    processing_results[database]["stdout"] = fo.database_directory + database + "_processing.out"
    processing_results[database]["stderr"] = fo.database_directory + tool + "_processing.err"

sys.stdout.write("Checking for tools required for the automated installation of third-party software.\n\n")
if subprocess.run(["which", args.pip], capture_output=True, env=env).returncode != 0 and not os.path.isfile(args.pip):
    pip_installed = False
    sys.stderr.write(args.pip + " not found. pip is required for installation of [snakemake, psutil, musicc, empanada]. These tools will not be installed during setup.\n")

if subprocess.run(["which", "make"], capture_output=True, env=env).returncode != 0:
    make_installed = False
    sys.stderr.write("make not found. Make is required for installation of [blast, bmtagger, markduplicates, diamond]. These tools will not be installed during setup.\n")

if subprocess.run(["which", "cmake"], capture_output=True, env=env).returncode != 0:
    cmake_installed = False
    sys.stderr.write("cmake not found. CMake is required for installation of [diamond]. These tools will not be installed during setup.\n")
sys.stdout.write("\n")

# Report any missing required tools
sys.stdout.write("Initial check for tools required for automated third-party software installation complete.\n\n")

# Try to install third party tools
sys.stdout.write("Beginning automated installation of third-party software.\n\n")
pip_install_list = []
if not args.no_snakemake:
    pip_install_list.append("snakemake")
    pip_install_list.append("psutil")
if not args.no_musicc:
    pip_install_list.append("musicc")
if not args.no_empanada:
    pip_install_list.append("empanada")
if len(pip_install_list) > 0:
    if pip_installed:
        sys.stdout.write("Installing python packages.\n\n")
        try:
            subprocess.run([args.pip, "install", "-t", fo.python_package_directory] + pip_install_list,
                           stdout=open(install_results["python_package"]["stdout"], "w"),
                           stderr=open(install_results["python_package"]["stderr"], "w"),
                           env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["python_package"]["error"] = True
            sys.stderr.write("Error: Failed to install: %s. Please see %s and %s for installation logs.\n" % (", ".join(pip_install_list), install_results["python_package"]["stdout"], install_results["python_package"]["stderr"]))

        # Adding additional checks because pip install failure doesn't exit with an error code
        if not args.no_snakemake and subprocess.run(["which", op.snakemake], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.snakemake):
            sys.stderr.write("Error: There was a problem installing snakemake locally via pip. Please see %s and %s for installation logs.\n" % (install_results["python_package"]["stdout"], install_results["python_package"]["stderr"]))

        if not args.no_musicc and subprocess.run(["which", config.steps.ortholog_abundance_correction.required_programs["musicc"]], capture_output=True, env=env).returncode != 0 and not os.path.isfile(config.steps.ortholog_abundance_correction.required_programs["musicc"]):
            sys.stderr.write("Error: There was a problem installing MUSiCC locally via pip. Please see %s and %s for installation logs.\n" % (install_results["python_package"]["stdout"], install_results["python_package"]["stderr"]))

        if not args.no_empanada and subprocess.run(["which", config.steps.ortholog_aggregation.required_programs["empanada"]], capture_output=True, env=env).returncode != 0 and not os.path.isfile(config.steps.ortholog_aggregation.required_programs["empanada"]):
            sys.stderr.write("Error: There was a problem installing EMPANADA locally via pip. Please see %s and %s for installation logs.\n" % (install_results["python_package"]["stdout"], install_results["python_package"]["stderr"]))

    else:
        sys.stderr.write("Warning: Pip not available, skipping installation of: %s.\n" % ", ".join(pip_install_list))
sys.stdout.write("\n")


def blast_executables_exist():
    return subprocess.run(["which", config.steps.host_filter.required_programs["blastn"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.host_filter.required_programs["blastn"])


if not args.no_blast and not blast_executables_exist():
    if make_installed:
        source_exists = True
        if not os.path.isdir(fo.source_directory + "ncbi-blast-2.2.31+-src"):
            sys.stdout.write("Downloading BLAST+ source.\n")
            try:
                subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-src.tar.gz", "-P", fo.source_directory],
                               stdout=open(install_results["blast"]["stdout"], "w"),
                               stderr=open(install_results["blast"]["stderr"], "w"),
                               env=env)
                subprocess.run(["tar", "-zxf", "ncbi-blast-2.2.31+-src.tar.gz"],
                               cwd=fo.source_directory,
                               stdout=open(install_results["blast"]["stdout"], "w"),
                               stderr=open(install_results["blast"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                install_error = True
                install_results["blast"]["error"] = True
                sys.stderr.write("Error: There was a problem downloading and/or expanding the BLAST+ source in preparation for local installation. Please see %s and %s for installation logs.\n" % (install_results["blast"]["stdout"], install_results["blast"]["stderr"]))
                source_exists = False

        if source_exists:
            sys.stdout.write("Installing BLAST+ locally.\n")
            try:
                subprocess.run(["./configure"],
                               cwd=fo.source_directory + "ncbi-blast-2.2.31+-src/c++/",
                               stdout=open(install_results["blast"]["stdout"], "w"),
                               stderr=open(install_results["blast"]["stderr"], "w"),
                               env=env)
                subprocess.run(["make", "all_r"],
                               cwd=fo.source_directory + "ncbi-blast-2.2.31+-src/c++/ReleaseMT/build/",
                               stdout=open(install_results["blast"]["stdout"], "w"),
                               stderr=open(install_results["blast"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                install_error = True
                install_results["blast"]["error"] = True
                sys.stderr.write("Error: There was a problem installing BLAST+ locally from source. Please see %s and %s for installation logs.\n" % (install_results["blast"]["stdout"], install_results["blast"]["stderr"]))

            # Adding additional check because BLAST make failure doesn't exit with an error code
            if not blast_executables_exist():
                sys.stderr.write("Error: There was a problem installing BLAST+ locally from source. Please see %s and %s for installation logs.\n" % (install_results["blast"]["stdout"], install_results["blast"]["stderr"]))

    else:
        sys.stderr.write("Warning: Make not available, skipping installation of: BLAST+.\n")
sys.stdout.write("\n")


def bmtools_executables_exist():
    return (subprocess.run(["which", config.steps.host_filter.required_programs["bmtool"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.host_filter.required_programs["bmtool"])) and (subprocess.run(["which", config.steps.host_filter.required_programs["bmfilter"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.host_filter.required_programs["bmfilter"])) and (subprocess.run(["which", config.steps.host_filter.required_programs["bmtagger"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.host_filter.required_programs["bmtagger"])) and (subprocess.run(["which", config.steps.host_filter.required_programs["extract_fa"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.host_filter.required_programs["extract_fa"]))


def srprism_executables_exist():
    return subprocess.run(["which", config.steps.host_filter.required_programs["srprism"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.host_filter.required_programs["srprism"])


if not args.no_bmtagger and (not bmtools_executables_exist() or not srprism_executables_exist()):
    if make_installed:
        if not bmtools_executables_exist():
            source_exists = True
            if not os.path.isdir(fo.source_directory + "bmtools/"):
                sys.stdout.write("Downloading BMTagger source.\n")
                try:
                    subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/bmtools.tar.gz", "-P", fo.source_directory],
                                   stdout=open(install_results["bmtools"]["stdout"], "w"),
                                   stderr=open(install_results["bmtools"]["stderr"], "w"),
                                   env=env)
                    subprocess.run(["tar", "-zxf", "bmtools.tar.gz"],
                                   cwd=fo.source_directory,
                                   stdout=open(install_results["bmtools"]["stdout"], "w"),
                                   stderr=open(install_results["bmtools"]["stderr"], "w"),
                                   env=env)
                except (EnvironmentError, subprocess.CalledProcessError):
                    install_error = True
                    install_results["bmtools"]["error"] = True
                    sys.stderr.write("Error: There was a problem downloading and/or expanding the BMTagger source in preparation for local installation. Please see %s and %s for installation logs.\n" % (install_results["bmtools"]["stdout"], install_results["bmtools"]["stderr"]))
                    source_exists = False

            if source_exists:
                sys.stdout.write("Installing BMTagger locally.\n")
                try:
                    subprocess.run(["make"],
                                   cwd=fo.source_directory + "bmtools/",
                                   stdout=open(install_results["bmtools"]["stdout"], "w"),
                                   stderr=open(install_results["bmtools"]["stderr"], "w"),
                                   env=env)
                except (EnvironmentError, subprocess.CalledProcessError):
                    install_error = True
                    install_results["bmtools"]["error"] = True
                    sys.stderr.write("Error: There was a problem installing BMTagger locally from source. Please see %s and %s for installation logs.\n" % (install_results["bmtools"]["stdout"], install_results["bmtools"]["stderr"]))

                # Adding additional check just in case the installation doesn't exit with an error code
                if not bmtools_executables_exist():
                    sys.stderr.write("Error: There was a problem installing BMTagger locally from source. Please see %s and %s for installation logs.\n" % (install_results["bmtools"]["stdout"], install_results["bmtools"]["stderr"]))

        if not srprism_executables_exist():
            source_exists = True
            if not os.path.isdir(fo.source_directory + "srprism/"):
                sys.stdout.write("Downloading srprism source.\n")
                try:
                    subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/src/srprism.tar.gz", "-P", fo.source_directory],
                                   stdout=open(install_results["srprism"]["stdout"], "w"),
                                   stderr=open(install_results["srprism"]["stderr"], "w"),
                                   env=env)
                    subprocess.run(["tar", "-zxf", "srprism.tar.gz"],
                                   cwd=fo.source_directory,
                                   stdout=open(install_results["srprism"]["stdout"], "w"),
                                   stderr=open(install_results["srprism"]["stderr"], "w"),
                                   env=env)
                except (EnvironmentError, subprocess.CalledProcessError):
                    install_error = True
                    install_results["srprism"]["error"] = True
                    sys.stderr.write("Error: There was a problem downloading and/or expanding the srprism source in preparation for local installation. Please see %s and %s for installation logs.\n" % (install_results["srprism"]["stdout"], install_results["srprism"]["stderr"]))
                    source_exists = False

            if source_exists:
                sys.stdout.write("Installing srprism locally.\n")
                try:
                    subprocess.run(["./ac.sh"],
                                   cwd=fo.source_directory + "srprism/gnuac/",
                                   stdout=open(install_results["srprism"]["stdout"], "w"),
                                   stderr=open(install_results["srprism"]["stderr"], "w"),
                                   env=env)
                    subprocess.run(["./configure"],
                                   cwd=fo.source_directory + "srprism/gnuac/",
                                   stdout=open(install_results["srprism"]["stdout"], "w"),
                                   stderr=open(install_results["srprism"]["stderr"], "w"),
                                   env=env)
                    subprocess.run(["make"],
                                   cwd=fo.source_directory + "srprism/gnuac/",
                                   stdout=open(install_results["srprism"]["stdout"], "w"),
                                   stderr=open(install_results["srprism"]["stderr"], "w"),
                                   env=env)
                except (EnvironmentError, subprocess.CalledProcessError):
                    install_error = True
                    install_results["srprism"]["error"] = True
                    sys.stderr.write("Error: There was a problem installing srprism locally from source. Please see %s and %s for installation logs.\n" % (install_results["srprism"]["stdout"], install_results["srprism"]["stderr"]))

                # Adding additional check just in case the installation doesn't exit with an error code
                if not srprism_executables_exist():
                    sys.stderr.write("Error: There was a problem installing srprism locally from source. Please see %s and %s for installation logs.\n" % (install_results["srprism"]["stdout"], install_results["srprism"]["stderr"]))

    else:
        sys.stderr.write("Warning: Make not available, skipping installation of: BMTagger.\n")
sys.stdout.write("\n")


def picard_executables_exist():
    return os.path.isfile(fo.source_directory + "picard.jar")


def samtools_executables_exist():
    return subprocess.run(["which", config.steps.duplicate_filter.required_programs["samtools"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.duplicate_filter.required_programs["samtools"])


if not args.no_markduplicates and (not picard_executables_exist() or not samtools_executables_exist()):
    if not picard_executables_exist():
        sys.stdout.write("Downloading PICARD jar.\n")
        try:
            subprocess.run(["wget", "https://github.com/broadinstitute/picard/releases/download/2.20.2/picard.jar", "-P", fo.source_directory],
                           stdout=open(install_results["picard"]["stdout"], "w"),
                           stderr=open(install_results["picard"]["stderr"], "w"),
                           env=env)
        except (EnvironmentError, subprocess.CalledProcessError):
            install_error = True
            install_results["picard"]["error"] = True
            sys.stderr.write("Error: There was a problem downloading the PICARD jar. Please see %s and %s for installation logs.\n" % (install_results["picard"]["stdout"], install_results["picard"]["stderr"]))

        # Adding additional check just in case the installation doesn't exit with an error code
        if not picard_executables_exist():
            sys.stderr.write("Error: There was a problem downloading the PICARD jar. Please see %s and %s for installation logs.\n" % (install_results["picard"]["stdout"], install_results["picard"]["stderr"]))

    if not samtools_executables_exist():
        if make_installed:
            source_exists = True
            if not os.path.isdir(fo.source_directory + "samtools-1.9/"):
                sys.stdout.write("Downloading samtools source.\n")
                try:
                    subprocess.run(["wget", "https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2", "-P", fo.source_directory],
                                   stdout=open(install_results["samtools"]["stdout"], "w"),
                                   stderr=open(install_results["samtools"]["stderr"], "w"),
                                   env=env)
                    subprocess.run(["tar", "-jxf", "samtools-1.9.tar.bz2"],
                                   cwd=fo.source_directory,
                                   stdout=open(install_results["samtools"]["stdout"], "w"),
                                   stderr=open(install_results["samtools"]["stderr"], "w"),
                                   env=env)
                except (EnvironmentError, subprocess.CalledProcessError):
                    install_error = True
                    install_results["samtools"]["error"] = True
                    sys.stderr.write("Error: There was a problem downloading and/or expanding the samtools source in preparation for local installation. Please see %s and %s for installation logs.\n" % (install_results["samtools"]["stdout"], install_results["samtools"]["stderr"]))
                    source_exists = False

            if source_exists:
                sys.stdout.write("Installing samtools locally.\n")
                try:
                    subprocess.run(["./configure", "--without-curses", "--disable-lzma"],
                                   cwd=fo.source_directory + "samtools-1.9/",
                                   stdout=open(install_results["samtools"]["stdout"], "w"),
                                   stderr=open(install_results["samtools"]["stderr"], "w"),
                                   env=env)

                    subprocess.run(["make"],
                                   cwd=fo.source_directory + "samtools-1.9/",
                                   stdout=open(install_results["samtools"]["stdout"], "w"),
                                   stderr=open(install_results["samtools"]["stderr"], "w"),
                                   env=env)
                except (EnvironmentError, subprocess.CalledProcessError):
                    install_error = True
                    install_results["samtools"]["error"] = True
                    sys.stderr.write("Error: There was a problem installing samtools locally from source. Please see %s and %s for installation logs.\n" % (install_results["samtools"]["stdout"], install_results["samtools"]["stderr"]))

                # Adding additional check just in case the installation doesn't exit with an error code
                if not samtools_executables_exist():
                    sys.stderr.write("Error: There was a problem installing samtools locally from source. Please see %s and %s for installation logs.\n" % (install_results["samtools"]["stdout"], install_results["samtools"]["stderr"]))

        else:
            sys.stderr.write("Warning: Make not available, skipping installation of: samtools.\n")
sys.stdout.write("\n")


def trimmomatic_executables_exist():
    return subprocess.run(["which", config.steps.quality_filter.required_programs["trimmer"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.quality_filter.required_programs["trimmer"])


if not args.no_trimmomatic and not trimmomatic_executables_exist():
    sys.stdout.write("Downloading TRIMMOMATIC jar.\n")
    try:
        subprocess.run(["wget", "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip", "-P", fo.source_directory],
                       stdout=open(install_results["trimmomatic"]["stdout"], "w"),
                       stderr=open(install_results["trimmomatic"]["stderr"], "w"),
                       env=env)
        subprocess.run(["unzip", fo.source_directory + "Trimmomatic-0.39.zip", "-d", fo.source_directory],
                       stdout=open(install_results["trimmomatic"]["stdout"], "w"),
                       stderr=open(install_results["trimmomatic"]["stderr"], "w"),
                       env=env)
    except (EnvironmentError, subprocess.CalledProcessError):
        install_error = True
        install_error["trimmomatic"]["error"] = True
        sys.stderr.write("Error: There was a problem downloading the TRIMMOMATIC jar. Please see %s and %s for installation logs.\n" % (install_results["trimmomatic"]["stdout"], install_results["trimmomatic"]["stderr"]))

    # Adding additional check just in case the installation doesn't exit with an error code
    if not trimmomatic_executables_exist():
        sys.stderr.write("Error: There was a problem downloading the TRIMMOMATIC jar. Please see %s and %s for installation logs.\n" % (install_results["trimmomatic"]["stdout"], install_results["trimmomatic"]["stderr"]))
sys.stdout.write("\n")


def diamond_executables_exist():
    return subprocess.run(["which", config.steps.map_reads_to_genes.required_programs["diamond"]], capture_output=True, env=env).returncode == 0 or os.path.isfile(config.steps.map_reads_to_genes.required_programs["diamond"])


if not args.no_diamond and not diamond_executables_exist():
    if make_installed and cmake_installed:
        source_exists = True
        if not os.path.isdir(fo.source_directory + "diamond-0.9.22/"):
            sys.stdout.write("Downloading DIAMOND source.\n")
            try:
                subprocess.run(["wget", "http://github.com/bbuchfink/diamond/archive/v0.9.22.tar.gz", "-P", fo.source_directory],
                               stdout=open(install_results["diamond"]["stdout"], "w"),
                               stderr=open(install_results["diamond"]["stderr"], "w"),
                               env=env)
                subprocess.run(["tar", "-zxf", "v0.9.22.tar.gz"],
                               cwd=fo.source_directory,
                               stdout=open(install_results["diamond"]["stdout"], "w"),
                               stderr=open(install_results["diamond"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                install_error = True
                install_results["diamond"]["error"] = True
                sys.stderr.write("Error: There was a problem downloading and/or expanding the DIAMOND source in preparation for local installation. Please see %s and %s for installation logs.\n" % (install_results["diamond"]["stdout"], install_results["diamond"]["stderr"]))
                source_exists = False

        if source_exists:
            sys.stdout.write("Installing DIAMOND locally.\n")
            try:
                subprocess.run(["cmake", ".", "-DCMAKE_INSTALL_PREFIX=" + os.getcwd() + "/" + fo.source_directory + "diamond-0.9.22/"],
                               cwd=fo.source_directory + "diamond-0.9.22/",
                               stdout=open(install_results["diamond"]["stdout"], "w"),
                               stderr=open(install_results["diamond"]["stderr"], "w"),
                               env=env)
                subprocess.run(["make", "install"],
                               cwd=fo.source_directory + "diamond-0.9.22/",
                               stdout=open(install_results["diamond"]["stdout"], "a"),
                               stderr=open(install_results["diamond"]["stderr"], "a"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                install_error = True
                install_results["diamond"]["error"] = True
                sys.stderr.write("Error: There was a problem installing DIAMOND locally from source. Please see %s and %s for installation logs.\n" % (install_results["diamond"]["stdout"], install_results["diamond"]["stderr"]))

            # Adding additional check just in case the installation doesn't exit with an error code
            if not diamond_executables_exist():
                sys.stderr.write("Error: There was a problem installing DIAMOND locally from source. Please see %s and %s for installation logs.\n" % (install_results["diamond"]["stdout"], install_results["diamond"]["stderr"]))

    else:
        sys.stderr.write("Warning: Make and/or CMake not available, skipping installation of: DIAMOND.\n")
sys.stdout.write("\n")

# Record any tools that are missing
all_missing_programs = set()
if subprocess.run(["which", op.python], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.python):
    all_missing_programs.add(op.python)

if subprocess.run(["which", op.java], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.java):
    all_missing_programs.add(op.java)

if subprocess.run(["which", op.snakemake], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.snakemake):
    all_missing_programs.add(op.snakemake)

makeblastdb_location = os.path.dirname(config.steps.host_filter.required_programs["blastn"]) + "/makeblastdb"
if subprocess.run(["which", makeblastdb_location], capture_output=True, env=env).returncode != 0 and not os.path.isfile(makeblastdb_location):
    all_missing_programs.add(makeblastdb_location)

for importer, modname, ispkg in pkgutil.iter_modules(config.steps.__path__):
    required_programs = importlib.import_module(".".join(["config.steps", modname])).required_programs.values()

    # Check if each file exists or is an executable in the environment
    step_missing_programs = []
    for required_program in required_programs:
        program_exists = os.path.isfile(required_program)
        executable_exists = subprocess.run(["which", required_program], capture_output=True, env=env).returncode == 0
        if not program_exists and not executable_exists:
            step_missing_programs.append(required_program)

    # If a program is missing, inform the user and add it to the set of missing programs
    if len(step_missing_programs) > 0:
        all_missing_programs = all_missing_programs.union(set(step_missing_programs))
sys.stdout.write("Automated installation of third-party software complete.\n\n")

sys.stdout.write("Beginning processing of reference data.\n\n")
# Try to create and process default database files
if not args.no_human_reference:
    if config.steps.host_filter.required_programs["bmtool"] not in all_missing_programs and config.steps.host_filter.required_programs["srprism"] not in all_missing_programs and makeblastdb_location not in all_missing_programs:
        source_exists = True
        if not os.path.isfile(fo.database_directory + op.host_database + ".fa"):
            sys.stdout.write("Downloading human reference.\n")
            try:
                subprocess.run(["wget", "ftp://ftp.ncbi.nlm.nih.gov/pub/agarwala/bmtagger/%s.fa" % op.host_database, "-P", fo.database_directory],
                               stdout=open(processing_results["human_reference"]["stdout"], "w"),
                               stderr=open(processing_results["human_reference"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["human_reference"]["error"] = True
                sys.stderr.write("Error: There was a problem downloading the human reference database in preparation for processing. Please see %s and %s for processing logs.\n" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"]))
                source_exists = False

        if not os.path.isfile(op.host_bitmask_file) and source_exists:
            sys.stdout.write("Creating human reference bitmask.\n")
            try:
                subprocess.run([config.steps.host_filter.required_programs["bmtool"], "-d", op.host_database_file + ".fa", "-o", op.host_bitmask_file, "-A", "0", "-w", "18"],
                               stdout=open(processing_results["human_reference"]["stdout"], "w"),
                               stderr=open(processing_results["human_reference"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["human_reference"]["error"] = True
                sys.stderr.write("Error: There was a problem creating a bitmask for the human reference database. Please see %s and %s for processing logs.\n" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"]))

        if not os.path.isfile(op.host_index_file + ".idx") and source_exists:
            sys.stdout.write("Creating human reference index.\n")
            try:
                subprocess.run([config.steps.host_filter.required_programs["srprism"], "mkindex", "-i", op.host_database_file + ".fa", "-o", op.host_index_file, "-M", "7168"],
                               stdout=open(processing_results["human_reference"]["stdout"], "w"),
                               stderr=open(processing_results["human_reference"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["human_reference"]["error"] = True
                sys.stderr.write("Error: There was a problem creating an index for the human reference database. Please see %s and %s for processing logs.\n" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"]))

        if not os.path.isfile(op.host_database_file + op.blast_db_suffix) and source_exists:
            sys.stdout.write("Creating human reference BLAST database.\n")
            try:
                subprocess.run([makeblastdb_location, "-in", op.host_database_file + ".fa", "-dbtype", "nucl", "-out", op.host_database_file],
                               stdout=open(processing_results["human_reference"]["stdout"], "w"),
                               stderr=open(processing_results["human_reference"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["human_reference"]["error"] = True
                sys.stderr.write("Error: There was a problem creating a BLAST database for the human reference database. Please see %s and %s for processing logs.\n" % (processing_results["human_reference"]["stdout"], processing_results["human_reference"]["stderr"]))
    else:
        human_reference_missing_programs = []
        if config.steps.host_filter.required_programs["bmtool"] in all_missing_programs:
            human_reference_missing_programs.append("bmtool")
        if config.steps.host_filter.required_programs["srprism"] in all_missing_programs:
            human_reference_missing_programs.append("srprism")
        if makeblastdb_location in all_missing_programs:
            human_reference_missing_programs.append("makeblastdb")
        sys.stderr.write("Warning: %s not available, skipping download and processing of human reference.\n" % ", ".join(human_reference_missing_programs))
sys.stdout.write("\n")

if not args.no_ko_mappings:
    sys.stdout.write("Downloading ortholog-to-grouping mappings.\n")
    for ortholog_to_grouping_mapping in op.ortholog_to_grouping_mappings:
        if not os.path.isfile(fo.ortholog_to_grouping_directory + ortholog_to_grouping_mapping + op.ortholog_to_grouping_suffix):
            try:
                subprocess.run(["wget", "https://github.com/borenstein-lab/fishtaco/raw/master/fishtaco/data/" + ortholog_to_grouping_mapping, "-O", fo.ortholog_to_grouping_directory + ortholog_to_grouping_mapping + op.ortholog_to_grouping_suffix],
                               stdout=open(processing_results["ortholog_mapping"]["stdout"], "w"),
                               stderr=open(processing_results["ortholog_mapping"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["ortholog_mapping"]["error"] = True
                sys.stderr.write("Error: There was a problem downloading the %s ortholog-to-grouping mapping file. Please see %s and %s for processing logs.\n" % (ortholog_to_grouping_mapping, processing_results["ortholog_mapping"]["stdout"], processing_results["ortholog_mapping"]["stderr"]))
sys.stdout.write("\n")

if args.uniprot:
    if config.steps.map_reads_to_genes.required_programs["diamond"] not in all_missing_programs and op.python not in all_missing_programs:
        source_exists = True
        if not os.path.isfile(fo.database_directory + op.target_database + ".fasta.gz"):
            sys.stdout.write("Downloading UniProt UniRef90 database.\n")
            try:
                subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/%s/%s.fasta.gz" % (op.target_database, op.target_database), "-P", fo.database_directory],
                               stdout=open(processing_results["uniprot"]["stdout"], "w"),
                               stderr=open(processing_results["uniprot"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem downloading the UniProt target database in preparation for processing. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
                source_exists = False

        if not os.path.isfile(op.target_database_file) and source_exists:
            sys.stdout.write("Creating UniRef90 DIAMOND database.\n")
            try:
                subprocess.run([config.steps.map_reads_to_genes.required_programs["diamond"], "makedb", "--in", fo.database_directory + op.target_database + ".fasta.gz", "-d", fo.database_directory + op.target_database],
                               stdout=open(processing_results["uniprot"]["stdout"], "w"),
                               stderr=open(processing_results["uniprot"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem creating a DIAMOND database for the  UniProt target database. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))

        if not os.path.isfile(op.gene_normalization_file) and source_exists:
            sys.stdout.write("Creating UniRef90 gene length table.\n")
            try:
                subprocess.run([op.python, fo.source_directory + "create_gene_length_table.py", fo.database_directory + op.target_database + ".fasta.gz", "--output", op.gene_normalization_file],
                               stdout=open(processing_results["uniprot"]["stdout"], "w"),
                               stderr=open(processing_results["uniprot"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem creating a gene length table for the  UniProt target database. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))

        source_exists = True
        if not os.path.isfile(fo.database_directory + "idmapping.dat.gz"):
            sys.stdout.write("Downloading UniRef90 gene-to-ortholog mapping.\n")
            try:
                subprocess.run(["wget", "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz", "-P", fo.database_directory],
                               stdout=open(processing_results["uniprot"]["stdout"], "w"),
                               stderr=open(processing_results["uniprot"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write(
                    "Error: There was a problem downloading the UniProt gene-to-ortholog mapping in preparation for processing. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
                source_exists = False

        if not os.path.isfile(op.gene_to_ortholog_file) and source_exists:
            sys.stdout.write("Formatting UniRef90 gene-to-ortholog mapping.\n")
            try:
                subprocess.run([op.python, fo.source_directory + "create_uniref_gene_to_ortholog.py", fo.database_directory + "idmapping.dat.gz", op.target_database, op.target_ortholog, "--output", op.gene_to_ortholog_file],
                               stdout=open(processing_results["uniprot"]["stdout"], "w"),
                               stderr=open(processing_results["uniprot"]["stderr"], "w"),
                               env=env)
            except (EnvironmentError, subprocess.CalledProcessError):
                processing_error = True
                processing_results["uniprot"]["error"] = True
                sys.stderr.write("Error: There was a problem formatting the gene-to-ortholog table for the  UniProt target database. Please see %s and %s for processing logs.\n" % (processing_results["uniprot"]["stdout"], processing_results["uniprot"]["stderr"]))
    else:
        uniprot_missing_programs = []
        if config.steps.map_reads_to_genes.required_programs["diamond"] in all_missing_programs:
            uniprot_missing_programs.append("DIAMOND")
        if op.python in all_missing_programs:
            uniprot_missing_programs.append("python")
        sys.stderr.write("Warning: %s not available, skipping download and processing of UniProt target database.\n" % ", ".join(uniprot_missing_programs))
sys.stdout.write("\n")

sys.stdout.write("*" * 80 + "\n")
sys.stdout.write("SETUP COMPLETE\n")
sys.stdout.write("Automated installation and reference data processing has finished. Please see above for details on software that could not be automatically installed and/or supporting data files that could not be generated.\n\n")

# Report any tools that are missing
if subprocess.run(["which", op.python], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.python):
    sys.stderr.write("Warning: %s is not present in your current environment. This will be required for running Python scripts. Please either install Python or change the 'python' variable in the config.operation.py submodule to point to a valid Python executable.\n" % op.python)

if subprocess.run(["which", op.java], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.java):
    sys.stderr.write("Warning: %s is not present in your current environment. This will be required for running Java jars. Please either install Java or change the 'java' variable in the config.operation.py submodule to point to a valid Java executable.\n" % op.java)

if subprocess.run(["which", op.snakemake], capture_output=True, env=env).returncode != 0 and not os.path.isfile(op.snakemake):
    sys.stderr.write("Warning: %s is not present in your current environment. This will be required for running MetaLAFFA. Please either install Snakemake or change the 'snakemake' variable in the config.operation.py submodule to point to a valid snakemake executable. If there was an issue installing Snakemake via PIP, you may consider using a Conda environment for your MetaLAFFA installation. See the README FAQ for more details.\n" % op.snakemake)

makeblastdb_location = os.path.dirname(config.steps.host_filter.required_programs["blastn"]) + "/makeblastdb"
if subprocess.run(["which", makeblastdb_location], capture_output=True, env=env).returncode != 0 and not os.path.isfile(makeblastdb_location):
    sys.stderr.write("Warning: %s is not present and executable in your current environment. This will be required for processing the default human reference database for host filtering. Check for BLAST+ installation errors if you tried to install BLAST+ with the setup script.\n" % makeblastdb_location)

for importer, modname, ispkg in pkgutil.iter_modules(config.steps.__path__):
    required_programs = importlib.import_module(".".join(["config.steps", modname])).required_programs.values()

    # Check if each file exists or is an executable in the environment
    step_missing_programs = []
    for required_program in required_programs:
        program_exists = os.path.isfile(required_program)
        executable_exists = subprocess.run(["which", required_program], capture_output=True, env=env).returncode == 0
        if not program_exists and not executable_exists:
            step_missing_programs.append(required_program)

    # If a program is missing, inform the user and add it to the set of missing programs
    if len(step_missing_programs) > 0:
        sys.stderr.write("Warning: program(s) missing for the %s pipeline step. The step may not run correctly without the program(s):\n" % modname)
        for missing_program in step_missing_programs:
            sys.stderr.write("%s\n" % missing_program)

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
            configured.write(line.format(python=op.python, MetaLAFFA_directory=os.getcwd(), PYTHONPATH="{{PYTHONPATH}}"))

sys.stdout.write("End of MetaLAFFA setup\n")
