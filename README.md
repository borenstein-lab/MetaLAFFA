Tutorial for running metaLAFFA
================

-   [Quick Start Summary](#quick-start-summary)
    -   [Downloading metaLAFFA](#downloading-metalaffa)
    -   [Configuring metaLAFFA](#configuring-metalaffa)
        -   [Using a custom functionally-annotated gene database](#using-a-custom-functionally-annotated-gene-database)
        -   [Configuring metaLAFFA to submit jobs to a cluster](#configuring-metalaffa-to-submit-jobs-to-a-cluster)
    -   [Trying out metaLAFFA](#trying-out-metalaffa)
-   [Full-length Tutorial](#full-length-tutorial)
    -   [Overview](#overview)
    -   [metaLAFFA step descriptions](#metalaffa-step-descriptions)
        -   [1. Host read filtering](#host-read-filtering)
        -   [2. Duplicate read filtering](#duplicate-read-filtering)
        -   [3. Quality trimming and filtering](#quality-trimming-and-filtering)
        -   [4. Mapping reads to gene IDs](#mapping-reads-to-gene-ids)
        -   [5. Filtering read mapping hits](#filtering-read-mapping-hits)
        -   [6. Gene counting](#gene-counting)
        -   [7. Ortholog counting](#ortholog-counting)
        -   [8. Correcting orthology group abundances](#correcting-orthology-group-abundances)
        -   [9. Aggregating orthology groups into higher-level functional descriptions](#aggregating-orthology-groups-into-higher-level-functional-descriptions)
        -   [10. Summarizing results of steps in metaLAFFA](#summarizing-results-of-steps-in-metalaffa)
    -   [Downloading metaLAFFA](#downloading-metalaffa-1)
        -   [Installing third-party tools and processing default databases](#installing-third-party-tools-and-processing-default-databases)
    -   [Configuring metaLAFFA](#configuring-metalaffa-1)
        -   [config/file\_organization.py](#configfile_organization.py)
        -   [config/operation.py](#configoperation.py)
        -   [config/cluster.py](#configcluster.py)
        -   [config/steps/\*.py](#configsteps.py)
        -   [config/library\_functions.py](#configlibrary_functions.py)
        -   [Configuring metaLAFFA to use a custom annotated gene database](#configuring-metalaffa-to-use-a-custom-annotated-gene-database)
        -   [Configuring metaLAFFA to submit jobs to a cluster](#configuring-metalaffa-to-submit-jobs-to-a-cluster-1)
        -   [Using a different step in metaLAFFA as a starting point](#using-a-different-step-in-metalaffa-as-a-starting-point)
    -   [Running metaLAFFA](#running-metalaffa)
        -   [Prepare input data files](#prepare-input-data-files)
        -   [Starting metaLAFFA](#starting-metalaffa)
        -   [Locating final outputs](#locating-final-outputs)
        -   [Restarting metaLAFFA](#restarting-metalaffa)
-   [References](#references)

**Repository:** <https://github.com/engal/metaLAFFA>

**Publication:** TBD

Quick Start Summary
-------------------

Follow these steps to annotate two example gut microbiome samples using metaLAFFA.

1.  If not already installed, install [**Python3**](https://www.python.org/downloads/) (version 3.7 or greater is required).

2.  Download and set up the metaLAFFA.

3.  Try out metaLAFFA with your own data or download some publicly available data to play with.

### Downloading metaLAFFA

To download metaLAFFA, first navigate to the directory in which you would like to save the pipeline. You can then clone the git repository as follows:

    # Make a local clone of the metaLAFFA git repository.
    git clone git@github.com:engal/metaLAFFA.git

Start pipeline setup by running the setup.py script:

    # Perform pipeline setup.
    cd metaLAFFA
    python3 setup.py

**Note**: using default options, it can take a few hours to install all third-party tools and download and process all default databases, takes ~10G of RAM, and ~30G of free disk space.

### Configuring metaLAFFA

##### Using a custom functionally-annotated gene database

By default, the setup script does not download any target protein database to map reads against for annotation purposes. You can either enable the `--uniprot` flag for the setup script, which will download the UniRef90 database (takes ~80G of free disk space) and create the necessary gene-to-ortholog and gene count normalization files automatically, or you can use your own existing database by modifying the pipeline configuration module. Using your own database (starting from a FASTA file of gene sequences) with the default pipeline configuration will require the following steps:

1.  Change the configuration module to point to your database and use your ortholog nomenclature of choice by editing the `config/operation.py` submodule, changing line 77 to:

        target_database = <name_of_your_database>

    and changing line 82 to:

        target_ortholog = <name_of_your_ortholog>

2.  Using DIAMOND, create a DIAMOND database from your database FASTA file:

        # Create a DIAMOND database
        diamond makedb --in <path_to_your_database_fasta> -d databases/<name_of_your_database>

    This will create a file called `databases/<name_of_your_database>.dmnd` what will then be used by DIAMOND to map reads to genes.

3.  Create a gene length normalization table from your database FASTA file using the included `create_gene_length_table.py` script:

        python3 src/create_gene_length_table.py <path_to_your_database_fasta> --output gene_normalization_files/<name_of_your_database>.norm

4.  Create a mapping file linking genes to orthologs. Annotations for databases can be in many formats, which means we are unable to provide a general-use script for producing a gene-to-ortholog mapping table. Instead, you will need to create a tab-separated table in the following format:

        ortholog_id1   gene_id1
        ortholog_id2   gene_id2
        ortholog_id1   gene_id3
        ...            ...

    where the first column consists of ortholog IDs and the second column consists of the IDs (same as in the FASTA) for genes associated with those orthologs. Note that this mapping can be many-to-many, with one ortholog being associated with multiple genes and one gene being associated with multiple orthologs. To work with default pipeline settings, place this gene-to-ortholog mapping table in `gene_to_ortholog_maps/<name_of_your_database>_to_<name_of_your_ortholog>.map`.

5.  If you are using KOs as your ortholog nomenclature, then you do not need to make any further changes for the pipeline to run and can ignore this step. However, if you are using a different functional annotation system, then you will need to make the following additional modifications.

    1.  Change the configuration module to use relative abundance normalization for ortholog counts, rather than MUSiCC because MUSiCC only works with KO counts. You can do this by editing the `config/steps/ortholog_abundance_correction.py` submodule, changing line 50 to:

                "method": "relative",

    2.  Configure the pipeline to skip aggregating orthologs to higher-level functional descriptions by commenting out lines 48-50 in `pipeline_steps.txt`. Do this by adding a `#` at the beginning of each line.

##### Configuring metaLAFFA to submit jobs to a cluster

By default, metaLAFFA is able to interface with Sun Grid Engine (SGE) and HTCondor clusters. This is achieved via the use of Python job submission wrapper scripts, included in the `src/` directory (`src/sge_submission_wrapper.py` and `src/condor_submission_wrapper.py` respectively). If your cluster uses a different cluster management system, then you will need to create your own job submission wrapper by following these steps:

1.  Copy the appropriate example job submission wrapper script to serve as a template for your new wrapper.

        cp src/<sge|condor>_submission_wrapper.py src/<name_of_your_cluster_system>_submission_wrapper.py

    Which example script you should use as a template depends on how you parameterize jobs when submitting them to the cluster. If you request cluster resources (memory, cores, etc.) via command-line arguments and just provide the name of a shell script to run (i.e. SGE uses `qsub <name_of_shell_script>`), then you should use the SGE wrapper as a template. If you instead specify cluster resources and which script to run via a config file (i.e. HTCondor uses `condor_submit <name_of_config_file>`), then you should use the HTCondor wrapper as a template.

2.  Edit the template script to submit a job on your cluster using appropriate resource requests and properly providing the script to run. The beginning of the template script covers necessary library imports and job parameter parsing. The second half of the template script handles cluster-specific processing and job submission. You should only need to edit the second half of the template script to let metaLAFFA interface with your cluster. Here is a brief explanation of what the example scripts are doing, to help you understand what changes you may need to make to interface with your cluster:

    1.  Starting from the cluster-specific section of the SGE template wrapper script:

        1.  First, the script intializes the list of components that make up the command to submit the job to the cluster. On SGE, the basic submission command is `qsub`, but you'll want to change this to the correct submission command for your cluster. The second component of the base submission command is the name of the script to run for the cluster job.

        2.  Next, there is some example processing for multi-core jobs. For the cluster where this was developed, multi-core jobs must request memory-per-core, rather than total memory for the job, so the script must calculate this. At the end of these calculations, the `qsub` command-line option for requesting multiple cores is added to the base submission command.

        3.  Once multi-core specific processing has been handled, the script adds command-line options for requesting memory, running time, and which directory to run the script in.

        4.  After the common resource requests and job settings have been added, the script next checks whether it should include a request for local disk space on a cluster node and whether it should request that resources be set aside until the job can run. The former is important if you have limited disk space on your shared file system and need to process intermediate files locally on cluster nodes during pipeline steps. The latter is useful for ensuring the jobs with large resource requests get to run.

        5.  Finally, the script runs the final submission command on the command-line.

    2.  Starting from the cluster-specific section of the HTCondor wrapper script:

        1.  First, the script converts memory requests into the appropriate units for HTCondor memory requests.

        2.  Next, the script opens up a config file (in `submission_files/`) specific to the pipeline operation it is submitting to the cluster. To keep config files distinct between parallel jobs, the script names this config file based on the name of the pipeline step and the instance of the pipeline step it is running (i.e. which input file is being run through this pipeline step).

        3.  The script then populates the config file with resource requests and settings. These include, in order of their addition in the script: the location of the script to run, the number of cores to request for the job, and the amount of memory to reserve for the job. The added `Queue` terminates the configuration file.

        4.  Finally, the script runs `condor_submit` on the command-line, providing the path to the config file for this job.

### Trying out metaLAFFA

1.  Obtain metagenomic shotgun sequencing data in FASTQ format to annotate. This can either be your own data, or you can download publicly available data to try annotating (e.g. from [the Human Microbiome Project](https://portal.hmpdacc.org/)).

2.  Format the names of your intial FASTQs to follow the default naming convention recognized by metaLAFFA:

    Forward-read FASTQs should be named `<sample>.R1.fastq.gz`

    Reverse-read FASTQs should be named `<sample>.R2.fastq.gz`

    Singleton/unpaired FASTQs should be named `<sample>.S.fastq.gz`

3.  Place the FASTQs in the `data/` directory.

4.  Start metaLAFFA with the following command (some tips on running metaLAFFA below):

        snakemake -c "python3 src/<name_of_your_cluster_system>_submission_wrapper.py" -j 50 --latency-wait 60

    1.  Since Snakemake needs to run for the entire annotation process, we recommend that you begin a [**screen session**](https://ss64.com/bash/screen.html) to start the pipeline and then detach the screen so you do not need to keep your terminal open.

    2.  The `-j 50` option specifies how many jobs should be in your cluster's queue at any one time, which limits the number of jobs that will be running at any one time but also avoids flooding the job queue with the potentially enormous number of jobs that may be generated. You should modify this setting according to your cluster resources and job queue etiquette.

    3.  The `--latency-wait 60` option specifies that Snakemake should wait up to a minute after a job finishes running before checking for the existence of the expected output file. Snakemake does this to determine whether the job completed successfully, and clusters can often have some delay between jobs technically finishing and output files being detectable when shared filesystems are involved. You should modify this setting according to the expected latency for your cluster environment.

    4.  metaLAFFA, by default, requests 10G of RAM per job, along with some step-specific memory requests according to the needs of the various tools used in the pipeline. You can modify the default resource requests by editing lines 9-13 of `config/cluster.py`, and you can modify step-specific resource requests by editing the `cluster_params` Python dictionary of each step's configuration submodule in `config/steps/<name_of_pipelin_step>.py`.

Full-length Tutorial
--------------------

### Overview

metaLAFFA is a pipeline for annotating shotgun metagenomic data with abundances of functional orthology groups. This process consists of several steps to go from raw FASTQs (with sequencing adapters removed) to functional profiles:

1.  Host read filtering (e.g. removing human DNA)
2.  Duplicate read filtering
3.  Quality trimming and filtering
4.  Mapping reads to genes
5.  Mapping genes to functional orthology groups (e.g. KOs)
6.  Aggregating orthologs into higher-level groupings (e.g. pathways)

This tutorial includes:

-   A description of what occurs during each step of the metaLAFFA pipeline
-   Details on downloading metaLAFFA
-   Instructions for configuring metaLAFFA to run on your cluster
-   A walkthrough of how to annotate metagenomic shotgun samples (example data or your own) using metaLAFFA

### metaLAFFA step descriptions

##### 1. Host read filtering

The first step in metaLAFFA is host read filtering. Host reads should be removed to ensure that later analyses are only considering the functional potential of the microbial community of interest. This is particularly important in cases where samples might consist of a majority of host reads, as in human fecal samples from hosts with certain gastrointestinal disorders that cause excessive shedding of epithelial cells (e.g. ulcerative colitis). Host read filtering is usually performed by aligning reads to one or more reference genomes for the host species and removing reads with match qualities above a certain threshold. Based on the original Human Microbiome Project (Huttenhower et al. 2012) study protocols, metaLAFFA defaults to using BMTagger (Rotmistrovsky and Agarwala 2011) for host read filtering. The default host database is the human reference hs37.

##### 2. Duplicate read filtering

The second step in metaLAFFA is duplicate read filtering. Duplicate read filtering aims to remove reads that erroneously appear in a sample multiple times due to artifacts of sequencing (PCR duplicates or optical duplicates). An explanation of the cause for PCR duplicates can be found [**here**](http://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/). However, PCR duplicate filtering may not always be appropriate for metagenomic studies because the (correct and biologically meaningful) presence of highly similar/equivalent reads is possible in metagenomics, especially when a given species or strain is at high abundance (Nayfach and Pollard 2016). In these scenarios, there are in fact more copies of the same sequence present, and filtering out duplicate reads may dampen the signal of this higher copy number within the community. Optical duplicates instead arise when, after amplification on a flowcell, a large cluster of duplicated reads registers as two separate clusters, and these are safe to remove. The default duplicate filtering program is MarkDuplicates from the Picard set of bioinformatics tools (“Picard,” n.d.), and by default only optical duplicates are removed.

##### 3. Quality trimming and filtering

The final pre-processing step in metaLAFFA is quality trimming and filtering. The aim of quality trimming and filtering is to reduce sequencing noise by trimming bases with quality below a certain threshold from one or both ends of a read, as well as removing any reads that fall below a specific length threshold. For this step, metaLAFFA uses Trimmomatic (Bolger, Lohse, and Usadel 2014) by default, trimming reads based on Trimmomatic's maximum information criterion and filtering out trimmed reads shorter than 60bp.

##### 4. Mapping reads to gene IDs

Given filtered and quality-controlled FASTQs, metaLAFFA next maps these reads against a database of protein sequences corresponding to identified genes. Note that this (read-based annotation) method is one of two main annotation approaches, the other being assembly-based annotation, such as performed by MOCAT2 (Kultima et al. 2016). Read-based annotation works by assigning functional annotations to individual reads rather than by assembling contigs from reads, assigning functional annotations to ORFs identified in the contigs, and then mapping reads to the ORFs to get counts.

Read annotations are most often performed by aligning reads to a database of gene sequences that have prior functional associations. These alignments are often performed as translated alignments, i.e. doing a 6-frame translation of a read's nucleotide sequence to possible amino acid sequences and then aligning those amino acid sequences to the protein sequences associated with genes. metaLAFFA uses DIAMOND (Buchfink, Xie, and Huson 2015) as its default aligner due to its speed and built-in parallelization. The UniRef90 reference database of gene sequences (from UniProt (Consortium 2018)) is used by default, though other databases of annotated gene sequences could be used if available (e.g. KEGG).

##### 5. Filtering read mapping hits

After mapping reads to gene sequences, metaLAFFA filters these hits to improve annotation quality. There are several ways to filter hits that depend on how many best hits are kept per read and whether hits to genes without functional annotations are kept. The differences in these methods are discussed in Carr et al. (Carr and Borenstein 2014), which analyzes the trade-offs in sensitivity and specificity that come from different choices in hit filtering. By default, metaLAFFA uses a custom Python script to filter hits and uses the highest specificity hit filtering method, keeping only hits with the best e-value for each read and allowing best hits to include genes without functional annotations.

##### 6. Gene counting

Once the hits are filtered, metaLAFFA counts gene abundances based on the number of reads that map to each gene in the sequence database. There are two important considerations when counting gene abundances based on mapped reads. The first is count normalization based on gene length, i.e. longer genes are more likely to be sequenced than shorter genes, and so a longer gene will be more likely to have more reads map to it than a shorter gene even though both genes might be at the same abundance and so the counts of reads mapping to a gene should be normalized by the length of the gene. The second consideration is how to count genes when a read maps to more than one gene. This can be done in one of two ways, either by giving a count of 1 to each gene a read maps to or by having each read contribute a count of 1 in total but dividing a read's count among the genes a read maps to. The former method can lead to double-counting issues and so metaLAFFA uses a fractional gene counting method by default. Gene length normalization and fractional counting are performed by a custom Python script in default metaLAFFA operation.

##### 7. Ortholog counting

The next step after calculating gene abundances is to convert those gene abundances to abundances of annotated functions. Usually, these come in the form of functional orthology group abundances as defined by various functional annotation systmes (e.g. KEGG (Kanehisa et al. 2018), MetaCyc (Caspi et al. 2007), etc.). These functional orthology groups usually correspond to the function a protein corresponding to a single gene may perform. Similar to the gene counting step, some genes may map to multiple orthology groups in a given system, and so each gene can either contribute its whole abundance to each associated orthology group or fractional contributions of its abundance to each associated orthology group. Again, due to double-counting issues, by default metaLAFFA uses a fraction contribution method when counting ortholog abundances based on gene abundances. Also, metaLAFFA calculates abundances of KEGG orthology groups (KOs) from UniProt gene functional annotations by default. By default, orthology group abundance calculations are performed using a custom Python script.

##### 8. Correcting orthology group abundances

Once initial orthology group abundances have been calculated, it is usually necessary to correct them to allow for comparisons between samples. Specifically, orthology group abundances based on read counts from standard shotgun metagenomic seuqencing are inherently compositional in nature (i.e. read count-based abundances reflect relative, rather than absolute, orthology group abundances). To account for this issue and convert orthology group abundances into a form that allows for valid standard statistical comparisons between samples, metaLAFFA uses MUSiCC (Manor and Borenstein 2015) during default operation. MUSiCC uses genes identified as universally single-copy across bacterial genomes to transform read-based relative orthology group abundances into average orthology group copy number across genomes within a sample. By default, metaLAFFA uses MUSiCC's intra-sample correction with learned models for incorporating universal single-copy gene abundances into the correction process.

##### 9. Aggregating orthology groups into higher-level functional descriptions

The final step in metaLAFFA is to aggregate functional orthology group abundances into abundances at higher-level functional descriptions (e.g. modules or pathways). This can help with interpretability of results, since initial orthology group counting can often result in thousands of unique functions. One way to handle this issue is to aggregate initial functional orthology groups into higher-level functional classifications based on their functional similarity or shared membership in a biological process (e.g. carbohydrate metabolism). This can result in a more manageable number of functions to analyze and can also reveal interesting trends in higher-level functional categories. By default, metaLAFFA uses a similar counting scheme to gene and orthology group counting when aggregating orthology groups, assigning an orthology group's abundance fractionally among all of the functional categories it belongs to at a given functional resolution. Associations between KOs and higher-level functional descriptions (modules and pathways) are determined from the KEGG BRITE hierarchy (Kanehisa et al. 2018) under default settings. During default operation, this is performed using a custom Python script.

##### 10. Summarizing results of steps in metaLAFFA

Another important component of metaLAFFA is the automatic summary statistic generation that accompanies each step in the pipeline. Default summary statistics include number of reads, average read length, and average base quality for input FASTQs and the resulting FASTQs after each FASTQ filtering step, as well as number of hits, number of mapped reads, number of associated genes, and number of annotated orthology groups as determined after each step in the pipeline. These summary statistics are merged at the end into a single table summarizing the results of each step of the pipeline.

### Downloading metaLAFFA

metaLAFFA is a Snakemake pipeline configured via Python scripts that manages the functional annotation of shotgun metagenomic data via various Python scripts and third-party tools.

First, if you do not have Python3 (version 3.7 or greater) installed, download the installation executable from [**here**](https://www.python.org/downloads/) and run it.

To download metaLAFFA, first navigate to the directory in which you would like to save the pipeline. You can then clone the git repository as follows:

    # Make a local clone of the metaLAFFA git repository.
    git clone git@github.com:engal/metaLAFFA.git

This contains the Snakemake file defining the pipeline, various Python scripts to perform pipeline steps that are not handled by third-party software, a setup script to help configure the environment that metaLAFFA will run in, and a Python module that configures how the steps of the pipeline are run. This does not include various third-party tools or default databases used by the pipeline.

##### Installing third-party tools and processing default databases

By default, metaLAFFA is built to use several third-party tools and one or two default databases for various pipeline steps. The third-party tools include Snakemake, BLAST+, BMTagger, MarkDuplicates, Samtools, Trimmomatic, DIAMOND, and MUSiCC, and by default metaLAFFA downloads the hs37 human reference for host filtering. You can complete pipeline setup, (installing these tools and database locally within the metaLAFFA directory) by running the setup.py script:

    # Perform default pipeline setup.
    cd metaLAFFA
    python3 setup.py

The `setup.py` script has various command-line options to modify the default setup procedure. These include options to skip the installation of third-party tools that you may already have installed or skip downloading the default human reference if you have a different host reference database to use. These options are of the form `--no_<name_of_tool>` and you can view the full list of these options by using the `-h` option.

**Note**: by default, metaLAFFA does not download a target database of functionally annotated genes to map reads to, instead assuming that you will configure metaLAFFA to use your own database. However, metaLAFFA is configured by default to work with KO annotations for the UniRef90 database, and you can have metaLAFFA download the UniRef90 database from UniProt by using the `--uniprot` option to the `setup.py` script. However, this database requires ~80G and can take a long time to download.

### Configuring metaLAFFA

metaLAFFA is configured via a Python module partitioned into various submodules targeted at specific aspects of the pipeline and defined in the "config" directory:

##### config/file\_organization.py

This submodule defines the file structure of metaLAFFA. Some parameters that you may want to change include:

`initial_data_directory`: Where to find the input metagenomic shotgun FASTQs to annotate.

`output_directory`: Where to generate the output files for each step in the pipeline.

`final_output_directory`: Where to write the final, desired output files once the pipeline is finished.

`tmp_dir`: Where to write the intermediate files that are generated during some steps of the pipeline but are not kept after the step is done running. By default this is `/tmp/`, which will use local storage on a cluster node, because this will be faster for file reading and writing. However, these intermediate files can be very large (they usually include uncompressed versions of the eventual output for a step) and if you have limited local storage on each cluster node, you may want to change this to process intermediate files in the metaLAFFA directory by setting it tot `./`.

##### config/operation.py

This submodule defines global settings that metaLAFFA uses during pipeline operation. These include whether output files should be compressed, the naming convention for the different types of paired and unpaired FASTQ files, the host database to use for filtering, the target database to use for mapping reads to genes, the file to use for mapping genes to orthologs, the file to use for aggregating orthologs to higher-level functional groupings, and the locations of executables for programs that are used by multiple steps in the pipeline (i.e. Python and Java executables). Some parameters that you may want to change include:

`sample_list`: The location of a list of sample IDs that tell metaLAFFA which samples to annotate. By default, metaLAFFA will annotate all samples it finds in the initial data directory whose file names match the expected pattern (e.g. <sample_id>.R1.fastq.gz for a forward read FASTQ). Providing a list of samples will limit metaLAFFA to only annotate those samples. This may be useful if you want to prioritiate the annotation of specific samples and annotate them first or just annotate samples in batches generally.

`fastq_types`: The mapping between the types of FASTQs (i.e. forward paired-reads, reverse paired-reads, and singleton reads) and the expected file name pattern. metaLAFFA will look for samples in the initial data directory that follow the pattern `<sample_id>.<type_id>.fastq.gz`, and you can use this configuration setting to alter the `type_id` that metaLAFFA will accept.

`host_database`: The name of the database to use for host filtering. By default, metaLAFFA uses the hs37 human reference genome for removing human reads from metagenomic samples, but you may want to use a different host reference for filtering depending on your project. This name should be the part that comes before the `.n*` suffix for the BLAST database (i.e. one of the default host reference database files is `hs37.nsq`). *Note*: BMTagger also requires a bitmask and an index, generated by `bmtools` and `srprism` respectively, that are associated with the host reference database. If you do change the host reference database, you will also need to generate `bitmasks/<your_database_name>.bitmask` and `indices/<your_database_name>.srprism.*` as described by [the HMP instructions for BMTagger](https://www.hmpdacc.org/hmp/doc/HumanSequenceRemoval_SOP.pdf).

`target_database`: The name of the database to use for mapping reads to genes with functional annotations. By default, metaLAFFA is configured to use the UniRef90 database, but you may want to use a different target database. See below for detailed instructions on all steps required to configure metaLAFFA to use a custom target database.

##### config/cluster.py

This submodule defines the default resource requests and settings for cluster jobs. These include memory, maximum running time, number of cores, and whether cluster resources should be reserved until a job is able to run. You can configure these settings to better fit your cluster by modifying the `default_cluster_params` dictionary.

##### config/steps/\*.py

Each step in the pipeline is configured by the submodule of the same name. These step-specific submodules define everything about a step in the pipeline, including expected naming patterns for input files, the output file naming patterns, custom cluster resource requests, the parameters for the program that runs during the step (e.g. method for filtering hits after mapping reads to genes or method for correcting ortholog abundances), and, most importantly, the actual operations to perform during the step.

A step's operations are defined using a Python function, which takes as its arguments the input files to the step, the output files to generate, and any additional information as determined by Snakemake. The function then runs the command-line operations that make up the step via the `subprocess` Python module. These submodules have been written such that advanced users can implement their own functions to redefine pipeline steps and easily swap that function in to run instead of the default function.

Some important default step-specific cluster resource requests (set with 10's to 100's of millions of reads in each FASTQ in mind) include:

`duplicate_filter`: 40G of RAM for MarkDuplicates

`quality_filter`: 40G of RAM for Trimmomatic

`map_reads_to_genes`: 220G of RAM and 22 cores for DIAMOND run with `--block_size 36` and `--index_chunks 1`. If you need to reduce the memory request for DIAMOND, you will also need to adjust the `block_size` and `index_chunks` settings in the read mapping configuration submodule as these control how much of the reference database DIAMOND loads into memory at a time. For more information on adjusting these DIAMOND parameters, see the DIAMOND manual [**here**](https://github.com/bbuchfink/diamond/raw/master/diamond_manual.pdf).

`hit_filter`: 200G of local disk space for creating a filtered version of DIAMOND output and then gzipping it.

`hit_filter_combine`: 250G of local disk space for concatenating hit-filtered DIAMOND output from FASTQs belong to the same sample (i.e. concatenating the results from the forward paired-read, reverse paired-read, and singlton read FASTQs from the same sample) and then gzipping the output.

##### config/library\_functions.py

This configuration submodule contains functions used in multiple places elsewhere in metaLAFFA to standardize their definitions.

##### Configuring metaLAFFA to use a custom annotated gene database

By default, the setup script does not download any target protein database to map reads against for annotation purposes. You can use your own existing database by modifying the pipeline configuration module. Using your own database (starting from a FASTA file of gene sequences) with the default pipeline configuration will require the following steps:

1.  Change the configuration module to point to your database and use your ortholog nomenclature of choice by editing the `config/operation.py` submodule, changing line 77 to:

        target_database = <name_of_your_database>

    and changing line 82 to:

        target_ortholog = <name_of_your_ortholog>

2.  Using DIAMOND, create a DIAMOND database from your database FASTA file:

        # Create a DIAMOND database
        diamond makedb --in <path_to_your_database_fasta> -d databases/<name_of_your_database>

    This will create a file called `databases/<name_of_your_database>.dmnd` what will then be used by DIAMOND to map reads to genes.

3.  Create a gene length normalization table from your database FASTA file using the included `create_gene_length_table.py` script:

        python3 src/create_gene_length_table.py <path_to_your_database_fasta> --output gene_normalization_files/<name_of_your_database>.norm

4.  Create a mapping file linking genes to orthologs. Annotations for databases can be in many formats, which means we are unable to provide a general-use script for producing a gene-to-ortholog mapping table. Instead, you will need to create a tab-separated table in the following format:

        ortholog_id1   gene_id1
        ortholog_id2   gene_id2
        ortholog_id1   gene_id3
        ...            ...

    where the first column consists of ortholog IDs and the second column consists of the IDs (same as in the FASTA) for genes associated with those orthologs. Note that this mapping can be many-to-many, with one ortholog being associated with multiple genes and one gene being associated with multiple orthologs. To work with default pipeline settings, place this gene-to-ortholog mapping table in `gene_to_ortholog_maps/<name_of_your_database>_to_<name_of_your_ortholog>.map`.

5.  If you are using KOs as your ortholog nomenclature, then you do not need to make any further changes for the pipeline to run and can ignore this step. However, if you are using a different functional annotation system, then you will need to make the following additional modifications.

    1.  Change the configuration module to use relative abundance normalization for ortholog counts, rather than MUSiCC because MUSiCC only works with KO counts. You can do this by editing the `config/steps/ortholog_abundance_correction.py` submodule, changing line 50 to:

                "method": "relative",

    2.  Configure the pipeline to either skip aggregating orthologs to higher-level functional descriptions or provide your own ortholog-to-functional-grouping mapping files.
        1.  To skip ortholog aggregation and just have the pipeline output corrected ortholog abundances as the final output, comment out lines 48-50 in `pipeline_steps.txt` by puting a `#` at the beginning of each line.

        2.  To provide your own mapping files, you must create tab-delimited tables in the form

                <name_of_your_ortholog>    group_id1    group_id2    ...
                ortholog_id1               1            0            ...
                ortholog_id2               1            1            ...
                ...                        ...          ...          ...

            where the first column contains the IDs for each ortholog and the remaining columns are labeled by an ID for a higher-level functional grouping. Each cell at row i and column j indicates whether ortholog i is part of functional group j, with 1 indicating that the ortholog is in the group and 0 indicating it is not in the group (For a full-scale example, the pipeline setup script downloads files for mapping KOs to modules and pathways, which are located in `ortholog_to_grouping_maps/`).

            1.  Write your mapping file to `ortholog_to_grouping_maps/<name_of_your_mapping>.map`.

            2.  Modify the configuration module to find those mapping files by editing the `config/operation.py` submodule, changine line 92 to:

                    ortholog_to_grouping_mappings = [<name_of_your_mapping>]

                *Note* You can use any number of mapping files (e.g. if you want to map to multiple levels of higher level function description such as modules and pathways). Just add each mapping file to the `ortholog_to_grouping_maps/` directory and include the name of each mapping in the list of mappings at line 92 in the `config/operation.py` submodule.

##### Configuring metaLAFFA to submit jobs to a cluster

By default, metaLAFFA is able to interface with Sun Grid Engine (SGE) and HTCondor clusters. This is achieved via the use of Python job submission wrapper scripts, included in the `src/` directory (`src/sge_submission_wrapper.py` and `src/condor_submission_wrapper.py` respectively). If your cluster uses a different cluster management system, then you will need to create your own job submission wrapper by following these steps:

1.  Copy the appropriate example job submission wrapper script to serve as a template for your new wrapper.

        cp src/<sge|condor>_submission_wrapper.py src/<name_of_your_cluster_system>_submission_wrapper.py

    Which example script you should use as a template depends on how you parameterize jobs when submitting them to the cluster. If you request cluster resources (memory, cores, etc.) via command-line arguments and just provide the name of a shell script to run (i.e. SGE uses `qsub <name_of_shell_script>`), then you should use the SGE wrapper as a template. If you instead specify cluster resources and which script to run via a config file (i.e. HTCondor uses `condor_submit <name_of_config_file>`), then you should use the HTCondor wrapper as a template.

2.  Edit the template script to submit a job on your cluster using appropriate resource requests and properly providing the script to run. The beginning of the template script covers necessary library imports and job parameter parsing. The second half of the template script handles cluster-specific processing and job submission. You should only need to edit the second half of the template script to let metaLAFFA interface with your cluster. Here is a brief explanation of what the example scripts are doing, to help you understand what changes you may need to make to interface with your cluster:

    1.  Starting from the cluster-specific section of the SGE template wrapper script:

        1.  First, the script intializes the list of components that make up the command to submit the job to the cluster. On SGE, the basic submission command is `qsub`, but you'll want to change this to the correct submission command for your cluster. The second component of the base submission command is the name of the script to run for the cluster job.

        2.  Next, there is some example processing for multi-core jobs. For the cluster where this was developed, multi-core jobs must request memory-per-core, rather than total memory for the job, so the script must calculate this. At the end of these calculations, the `qsub` command-line option for requesting multiple cores is added to the base submission command.

        3.  Once multi-core specific processing has been handled, the script adds command-line options for requesting memory, running time, and which directory to run the script in.

        4.  After the common resource requests and job settings have been added, the script next checks whether it should include a request for local disk space on a cluster node and whether it should request that resources be set aside until the job can run. The former is important if you have limited disk space on your shared file system and need to process intermediate files locally on cluster nodes during pipeline steps. The latter is useful for ensuring the jobs with large resource requests get to run.

        5.  Finally, the script runs the final submission command on the command-line.

    2.  Starting from the cluster-specific section of the HTCondor wrapper script:

        1.  First, the script converts memory requests into the appropriate units for HTCondor memory requests.

        2.  Next, the script opens up a config file (in `submission_files/`) specific to the pipeline operation it is submitting to the cluster. To keep config files distinct between parallel jobs, the script names this config file based on the name of the pipeline step and the instance of the pipeline step it is running (i.e. which input file is being run through this pipeline step).

        3.  The script then populates the config file with resource requests and settings. These include, in order of their addition in the script: the location of the script to run, the number of cores to request for the job, and the amount of memory to reserve for the job. The added `Queue` terminates the configuration file.

        4.  Finally, the script runs `condor_submit` on the command-line, providing the path to the config file for this job.

##### Using a different step in metaLAFFA as a starting point

While metaLAFFA defaults to taking in unfiltered FASTQs, it is also possible for metaLAFFA to take input files appropriate for any step in the pipeline. For example, if you have already processed and filtered your FASTQ files, you can use those as input to the pipeline as long as they follow the correct naming scheme (e.g. `<sample_id>.<type_id>.fastq.gz` for FASTQs). To find the naming scheme for input files for different steps in the pipeline, you can look in the configuration submodule associated with the step and look at the `input` setting.

You must indicate which intermediate step you want metaLAFFA to start from by appropriately modifying the `pipeline_steps.txt` file. This file tells the pipeline what the input to each step in the pipeline is (either initial input data or the output of previous steps). For example, by default the host filtering step takes in the initial input data as its input, and this is indicated by the line `host_filter:INPUT` (step name on the left, followed by a colon and then a comma-separated list of where the input to the step comes from). If you want to start from FASTQs that you have already filtered, you would make the following changes to `pipeline_steps.txt`:

-   Change the reading mapping step input to be the initial input data (i.e. changing `map_reads_to_genes:quality_filter` -&gt;; `map_reads_to_genes:INPUT`).
-   Comment-out or delete steps that you don't want to run (i.e. for all lines for steps starting with `host_filter`, `duplicate_filter`, or `quality_filter`, either put a `#` at the front of the line or delete the line).

Once you've made those changes, metaLAFFA will start from the read mapping step and use files in the `data` directory as input. There are step markers that can be used to indicate certain features of each step, and these should be placed at the beginning of the line for the associated step in`pipelin_steps.txt`. These markers include:

`$`: This step produces a final output file that should be copied to a location apart from the rest of the output files to make it easier to locate. Additionally, to avoid unnecessary work, Snakemake will not run any steps that are not required to generate any final output files.

`*`: This step produces a summary table for one or more steps in the pipelin. At the end of the pipeline, the outputs of these steps are merged into a single master summary table as a final output.

### Running metaLAFFA

Here we provide a walkthrough for running metaLAFFA on a metagenomic shotgun sequencing dataset starting from unfiltered FASTQs.

##### Prepare input data files

1.  Obtain metagenomic shotgun sequencing data in FASTQ format to annotate. This can either be your own data, or you can download publicly available data to try annotating (e.g. from [the Human Microbiome Project](https://portal.hmpdacc.org/)).

2.  Format the names of your intial FASTQs to follow the default naming convention recognized by metaLAFFA:

    Forward-read FASTQs should be named `<sample>.R1.fastq.gz`

    Reverse-read FASTQs should be named `<sample>.R2.fastq.gz`

    Singleton/unpaired FASTQs should be named `<sample>.S.fastq.gz`

    **Note** metaLAFFA expects all three FASTQ types for each sample. If one or more do not exist for a sample (e.g. some samples only have paired reads while others also have some singletons), metaLAFFA will generate appropriately-named empty dummy files as placeholders for the missing expected files. The creation of these dummy files is logged in `dummy_input_files_generated.txt`. If you start metaLAFFA from a step that does not take FASTQs as input (e.g. you have your own DIAMOND output from a separate source and want to annotate it), metaLAFFA will do the same check for expected files and generate (and log the creation of) dummy files for any expected files that are missing.

3.  Place the FASTQs in the `data/` directory.

##### Starting metaLAFFA

1.  Start metaLAFFA with the following command (some tips on running metaLAFFA below):

        snakemake -c "python3 src/<name_of_your_cluster_system>_submission_wrapper.py" -j 50 --latency-wait 60

    1.  Since Snakemake needs to run for the entire annotation process, we recommend that you begin a [**screen session**](https://ss64.com/bash/screen.html) to start the pipeline and then detach the screen so you do not need to keep your terminal open.

    2.  The `-j 50` option specifies how many jobs should be in your cluster's queue at any one time, which limits the number of jobs that will be running at any one time but also avoids flooding the job queue with the potentially enormous number of jobs that may be generated. You should modify this setting according to your cluster resources and job queue etiquette.

    3.  The `--latency-wait 60` option specifies that Snakemake should wait up to a minute after a job finishes running before checking for the existence of the expected output file. Snakemake does this to determine whether the job completed successfully, and clusters can often have some delay between jobs technically finishing and output files being detectable when shared filesystems are involved. You should modify this setting according to the expected latency for your cluster environment.

    4.  metaLAFFA, by default, requests 10G of RAM per job, along with some step-specific memory requests according to the needs of the various tools used in the pipeline. You can modify the default resource requests by editing lines 9-13 of `config/cluster.py`, and you can modify step-specific resource requests by editing the `cluster_params` Python dictionary of each step's configuration submodule in `config/steps/<name_of_pipelin_step>.py`.

##### Locating final outputs

Once the pipeline has finished running, you can find your desired output files in the default final output directory, `output/final/`. If you ran the pipeline in the default configuration, these should include:

`output/final/KOvsMODULE_BACTERIAL_KEGG_2013_07_15.aggregated_orthologs.tab`: The module-level functional profile

`output/final/KOvsPATHWAY_BACTERIAL_KEGG_2013_07_15.aggregated_orthologs.tab`: The pathway-level functional profile

`output/final/orthologs.tab`: The average KO copy number per genome, as determined by MUSiCC

`output/final/summary.txt`: The master table summarizing each step in the pipeline

##### Restarting metaLAFFA

You can stop metaLAFFA at any point (e.g. CTRL-C while running Snakemake) and if you rerun the same Snakemake command, the pipeline will pick up where it left off.

References
----------

Bolger, Anthony M, Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A Flexible Trimmer for Illumina Sequence Data.” *Bioinformatics* 30 (15). Oxford University Press: 2114–20.

Buchfink, Benjamin, Chao Xie, and Daniel H Huson. 2015. “Fast and Sensitive Protein Alignment Using Diamond.” *Nature Methods* 12 (1). Nature Publishing Group: 59.

Carr, Rogan, and Elhanan Borenstein. 2014. “Comparative Analysis of Functional Metagenomic Annotation and the Mappability of Short Reads.” *PLoS One* 9 (8). Public Library of Science: e105776.

Caspi, Ron, Hartmut Foerster, Carol A Fulcher, Pallavi Kaipa, Markus Krummenacker, Mario Latendresse, Suzanne Paley, et al. 2007. “The Metacyc Database of Metabolic Pathways and Enzymes and the Biocyc Collection of Pathway/Genome Databases.” *Nucleic Acids Research* 36 (suppl\_1). Oxford University Press: D623–D631.

Consortium, UniProt. 2018. “UniProt: A Worldwide Hub of Protein Knowledge.” *Nucleic Acids Research* 47 (D1). Oxford University Press: D506–D515.

Huttenhower, Curtis, Dirk Gevers, Rob Knight, Sahar Abubucker, Jonathan H Badger, Asif T Chinwalla, Heather H Creasy, et al. 2012. “Structure, Function and Diversity of the Healthy Human Microbiome.” *Nature* 486 (7402). Nature Publishing Group: 207.

Kanehisa, Minoru, Yoko Sato, Miho Furumichi, Kanae Morishima, and Mao Tanabe. 2018. “New Approach for Understanding Genome Variations in Kegg.” *Nucleic Acids Research* 47 (D1). Oxford University Press: D590–D595.

Kultima, Jens Roat, Luis Pedro Coelho, Kristoffer Forslund, Jaime Huerta-Cepas, Simone S Li, Marja Driessen, Anita Yvonne Voigt, Georg Zeller, Shinichi Sunagawa, and Peer Bork. 2016. “MOCAT2: A Metagenomic Assembly, Annotation and Profiling Framework.” *Bioinformatics* 32 (16). Oxford University Press: 2520–3.

Manor, Ohad, and Elhanan Borenstein. 2015. “MUSiCC: A Marker Genes Based Framework for Metagenomic Normalization and Accurate Profiling of Gene Abundances in the Microbiome.” *Genome Biology* 16 (1). BioMed Central: 53.

Nayfach, Stephen, and Katherine S Pollard. 2016. “Toward Accurate and Quantitative Comparative Metagenomics.” *Cell* 166 (5). Elsevier: 1103–16.

“Picard.” n.d. <http://broadinstitute.github.io/picard/>.

Rotmistrovsky, Kirill, and Richa Agarwala. 2011. “BMTagger: Best Match Tagger for Removing Human Reads from Metagenomics Datasets.” *Unpublished*.
