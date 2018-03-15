# Snakemake pipeline for Borenstein lab metagenomics annotation pipeline
# Author: Adrian Verster
# December 2017

# snakemake -p -c "qsub {params.cluster}" -j 50 --latency-wait 60


  import config,gzip,os

if config.many_subdirectories:
    #If you want subdirectories
    join_char = "/"
else:
    #if you want long filenames
    join_char = "_"


#diamond_output_suffix = join_char.join([x for x in ["diamond_%s_sensitivity_%s" %( config.alignment_method, config.sensitivity), "top_percentage_%s" %( str(config.top_percentage)), "max_e_value_%s" %(str(config.max_e_value)) ] if x])
#hitfiltering_output_suffix = diamond_output_suffix + join_char + join_char.join([x for x in ["best_n_hits_%s" %( str(config.best_n_hits) ), "filtering_method_%s" %( config.filtering_method) ] if x])
#genemapper_output_suffix = hitfiltering_output_suffix + join_char + join_char.join([x for x in ["count_method_%s" %( config.count_method) ] if x])
#normalization_output_suffix = genemapper_output_suffix + join_char + join_char.join([x for x in ["norm_method_%s" %( config.norm_method) ,"musicc_correction_method_%s" %( config.musicc_correction_method)] if x])
#functionalsummary_output_suffix = normalization_output_suffix + join_char + join_char.join([x for x in ["mapping_matrix_%s" %( os.path.basename(config.mapping_matrix) ), "summary_method_%s" %( config.summary_method), "functional_level_%s" %(config.functional_level) ] if x])


#Acronyms for the parameters are as follows:
# D diamond alignment_method
# S diamond sensitivity
# TP diamond top_percentage
# MEV diamond max_e_value
# DB diamond database
# BNH hit filtering best_n_hits
# FM hit filtering filtering_method
# CMG gene mapper count_method_gene
# CMK ko mapper count_method_ko
# NM normalization norm_method
# MCM normalization musicc_correction_method
# MM functionalsummary mapping matrix
# SM functionalsummary summary method
# FL functionalsummary functional level

db_name = config.kegg_version + "-" + config.taxon
diamond_output_suffix = join_char.join([x for x in ["D_%s_S_%s" %( config.alignment_method, config.sensitivity), "TP_%s" %( str(config.top_percentage)), "MEV_%s" %(str(config.max_e_value)), "DB_%s" %(db_name) ] if x])
hitfiltering_output_suffix = diamond_output_suffix + join_char + join_char.join([x for x in ["BNH_%s" %( str(config.best_n_hits) ), "FM_%s" %( config.filtering_method) ] if x])
genemapper_output_suffix = hitfiltering_output_suffix + join_char + join_char.join([x for x in ["CMG_%s" %( config.count_method_gene) ] if x])
komapper_output_suffix = genemapper_output_suffix + join_char + join_char.join([x for x in ["CMK_%s" %(config.count_method_ko)] if x])
normalization_output_suffix = komapper_output_suffix + join_char + join_char.join([x for x in ["NM_%s" %( config.norm_method) ,"MCM_%s" %( config.musicc_correction_method)] if x])
functionalsummary_output_suffix = normalization_output_suffix + join_char + join_char.join([x for x in ["MM_%s" %( db_name ), "SM_%s" %( config.summary_method), "FL_%s" %(config.summary_level) ] if x])
default_cluster_params = "-cwd -l mfree=10G,h_rt=24:0:0" 
delete_intermediates = config.delete_intermediates

#Without this line snakemake will sometimes fail a job because it fails to detect the output file due to latency
#shell.suffix("; sleep 40")

if config.samples_oi is None:
    if len(os.listdir(config.fastq_directory)) > 0:
        SAMPLES = list(set([x.split(".")[0] for x in os.listdir(config.fastq_directory) if x.split(".")[0] != ""]))
    elif len(os.listdir(config.quality_filtered_directory)):
        SAMPLES = list(set([x.split(".")[0] for x in os.listdir(config.quality_filtered_directory) if x.split(".")[0] != ""]))
    else:
        sys.exit("No samples in either the config.fastq_directory or the config.quality_filtered_directory")
else:
    SAMPLES = []
    with open(config.samples_oi) as f:
        for line in f:
            line = line.rstrip("\n")
            SAMPLES.append(line)

#Only unique SAMPLES
SAMPLES = list(set(SAMPLES))

wildcard_constraints:
    sample="[A-Za-z0-9_]+",
    type="[A-Za-z0-9_]+"

rule all:
    input:
        #config.module_profiles_directory + "functionalsummary." + functionalsummary_output_suffix + ".gz",
        #config.summary_directory + "functional_summary_summary." + functionalsummary_output_suffix + ".txt",
        #config.module_profiles_directory + functionalsummary_output_suffix + "/functionalsummary.gz"
        #config.summary_directory +  genemapper_output_suffix + "/ko_mapper_summary.txt"
        #config.summary_directory + "duplicate_filter_summary.txt"
        #Summaries
        config.summary_directory + functionalsummary_output_suffix + "/Summaries_All.txt"
        #expand(config.summary_directory + "{sample}.fastq_summary.txt", sample = SAMPLES)
        #config.summary_directory + "fastq_summary.txt"


rule clean_all:
   run:
        shell( "rm -f -r %s4_blast_results/* && rm -f -r %s5_filtered_blast_results/* && rm -f -r %s6_gene_profiles/* && rm -f -r %s7_ko_profiles/* && rm -f -r %s9_module_profiles/* && rm -f -r %s8_normalized_ko_profiles/* && rm -f -r %s2_duplicate_filtered/* && rm -f -r %s1_host_filtered/* && rm -f -r %s3_quality_filtered/*" %(config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir ) )
        shell("rm -f -r %slogging/* && rm -f -r %ssummaries/*" %( config.output_dir,config.output_dir ) )

rule clean_sub:
    run:
        shell( "rm -f -r %s4_blast_results/* && rm -f -r %s5_filtered_blast_results/ && rm -f -r %s6_gene_profiles/* && rm -f -r %s7_ko_profiles/* && rm -f -r %s9_module_profiles/* && rm -f -r %s8_normalized_ko_profiles/* " %(config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir,config.output_dir ) )


rule fastq_summary:
    input:
        config.fastq_directory + "{sample}.{type}.fastq.gz",
    output:
        config.summary_directory + "{sample}.{type}.fastq_summary.txt",
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_input_fastq.py {input} > {output}"

def fastq_summary_combine1_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R1.fastq.gz".format(wildcards=wildcards)):
        out["R1"] = config.summary_directory + "{wildcards.sample}.R1.fastq_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R2.fastq.gz".format(wildcards=wildcards)):
        out["R2"] = config.summary_directory + "{wildcards.sample}.R2.fastq_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.S.fastq.gz".format(wildcards=wildcards)):
        out["S"] = config.summary_directory + "{wildcards.sample}.S.fastq_summary.txt".format(wildcards=wildcards)
    return out

rule fastq_summary_combine1:
    input:
        unpack(fastq_summary_combine1_DetermineFiles)
    output:
        config.summary_directory + "{sample}.fastq_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        shell("head -1 %s > %s" %(input[0], output ) )
        for i in range(len(input)):
            shell("tail -n +2 %s >> %s" %(input[i], output))
        shell("rm -f {input}")

rule fastq_summary_combine2:
    input:
        expand( config.summary_directory + "{sample}.fastq_summary.txt", sample = SAMPLES )
    output:
        config.summary_directory + "fastq_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)
        shell("rm -f {input}")


rule host_filter_paired:
    input:
        R1=config.fastq_directory + "{sample}.R1.fastq.gz",
        R2=config.fastq_directory + "{sample}.R2.fastq.gz"
    output:
        R1=config.host_filtered_directory + "{sample}.R1.hostfilter.fastq.gz",
        R2=config.host_filtered_directory + "{sample}.R2.hostfilter.fastq.gz",
        host_marked=config.host_filtered_directory + "{sample}.R.hostmarked.fastq.gz"
    params:
        cluster=default_cluster_params
    benchmark:
        "".join([config.log_directory, "host_filter_paired.{sample}.log"])
    run:
        out_F_nonzip = output.R1.rstrip(".gz")
        out_R_nonzip = output.R2.rstrip(".gz")
        out_host_nonzip = output.host_marked.rstrip(".gz")
        shell( "src/host_filtering_wrapper.sh {input.R1} %s %s --paired_fastq {input.R2} --paired_fastq_output %s" %(out_host_nonzip, out_F_nonzip, out_R_nonzip) )
        shell( "gzip %s" %(out_F_nonzip) )
        shell( "gzip %s" %(out_R_nonzip) )
        shell( "gzip %s" %(out_host_nonzip) )
        #Delete intermediate # Deleting the intermediates should also get rid of the marked read file
        if delete_intermediates:
            shell("rm -f {input.R1}")
            shell("rm -f {input.R2}")

rule host_filter_singleton:
    input:
        S=config.fastq_directory + "{sample}.S.fastq.gz",
    output:
        S=config.host_filtered_directory + "{sample}.S.hostfilter.fastq.gz",
        host_marked=config.host_filtered_directory + "{sample}.S.hostmarked.fastq.gz"
    params:
        cluster=default_cluster_params
    benchmark:
        "".join([config.log_directory, "host_filter_singleton.{sample}.log"])
    run:
        out_nonzip = output.S.rstrip(".gz")
        out_host_nonzip = output.host_marked.rstrip(".gz")
        shell(" src/host_filtering_wrapper.sh {input.S} %s %s" %(out_host_nonzip, out_nonzip) )
        shell( "gzip %s" %(out_nonzip) )
        shell( "gzip %s" %(out_host_nonzip) )
        #Delete intermediate # Deleting the intermediates should also get rid of the marked read file
        if delete_intermediates:
            shell("rm -f {input.S}")

rule host_filter_summary:
    input:
        config.host_filtered_directory + "{sample}.{type}.hostfilter.fastq.gz"
    output:
        config.summary_directory + "{sample}.{type}.hostfilter_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_host_filtering.py {input} > {output}"

def host_filter_summary_combine1_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R1.fastq.gz".format(wildcards=wildcards)):
        out["R1"] = config.summary_directory + "{wildcards.sample}.R1.hostfilter_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R2.fastq.gz".format(wildcards=wildcards)):
        out["R2"] = config.summary_directory + "{wildcards.sample}.R2.hostfilter_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.S.fastq.gz".format(wildcards=wildcards)):
        out["S"] = config.summary_directory + "{wildcards.sample}.S.hostfilter_summary.txt".format(wildcards=wildcards)
    return out

rule host_filter_summary_combine1:
    input:
        unpack(host_filter_summary_combine1_DetermineFiles)
    output:
        config.summary_directory + "{sample}.hostfilter_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        shell("head -1 %s > %s" %(input[0], output ) )
        for i in range(len(input)):
            shell("tail -n +2 %s >> %s" %(input[i], output))
        shell("rm -f {input}")

rule host_filter_summary_combine2:
    input:
        expand( config.summary_directory + "{sample}.hostfilter_summary.txt", sample = SAMPLES )
    output:
        config.summary_directory + "hostfilter_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)
        shell("rm -f {input}")

rule duplicate_filter_singleton:
    input:
        S=config.host_filtered_directory + "{sample}.S.hostfilter.fastq.gz",
    output:
        S=config.duplicate_filtered_directory + "{sample}.S.dupfilter.fastq.gz",
    params:
        cluster=default_cluster_params
    benchmark:
        "".join([config.log_directory, "duplicate_filter_singleton.{sample}.log"])
    run:
        out_nonzip = output.S.rstrip(".gz")
        shell(" src/duplicate_filtering_wrapper.sh {input.S} {wildcards.sample} dummy1 dummy2 %s " %(out_nonzip) )
        shell( "gzip %s" %(out_nonzip) )
        #Delete intermediate # Deleting the intermediates should also get rid of the marked read file
        if delete_intermediates:
            shell("rm -f {input.S}")

rule duplicate_filter_paired:
    input:
        R1=config.host_filtered_directory + "{sample}.R1.hostfilter.fastq.gz",
        R2=config.host_filtered_directory + "{sample}.R2.hostfilter.fastq.gz",
    output:
        R1=config.duplicate_filtered_directory + "{sample}.R1.dupfilter.fastq.gz",
        R2=config.duplicate_filtered_directory + "{sample}.R2.dupfilter.fastq.gz",
        metrics_output=config.duplicate_filtered_directory + "{sample}.R.metricout.fastq.gz",
        duplicate_marked=config.duplicate_filtered_directory + "{sample}.R.duplicatemarked.fastq.gz",
    params:
        cluster="-cwd -l mfree=32G,h_rt=24:0:0"
    benchmark:
        "".join([config.log_directory, "duplicate_filter_paired.{sample}.log"])
    run:
        out_F_nonzip = output.R1.rstrip(".gz")
        out_R_nonzip = output.R2.rstrip(".gz")
        out_metrics_nonzip = output.metrics_output.rstrip(".gz")
        out_duplicate_nonzip = output.duplicate_marked.rstrip(".gz")
        shell( "src/duplicate_filtering_wrapper.sh {input.R1} {wildcards.sample} %s %s %s --paired_fastq {input.R2} --paired_fastq_output %s" %(out_metrics_nonzip, out_duplicate_nonzip,out_F_nonzip, out_R_nonzip) )
        shell( "gzip %s" %(out_F_nonzip) )
        shell( "gzip %s" %(out_R_nonzip) )
        shell( "gzip %s" %(out_metrics_nonzip) )
        shell( "gzip %s" %(out_duplicate_nonzip) )
        #Delete intermediate # Deleting the intermediates should also get rid of the marked read file
        if delete_intermediates:
            shell("rm -f {input.R1}")
            shell("rm -f {input.R2}")

rule duplicate_filter_summary:
    input:
        config.duplicate_filtered_directory + "{sample}.{type}.dupfilter.fastq.gz"
    output:
        config.summary_directory + "{sample}.{type}.duplicate_filter_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_duplicate_filtering.py {input} > {output}"

def duplicate_filter_summary_combine1_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R1.fastq.gz".format(wildcards=wildcards)):
        out["R1"] = config.summary_directory + "{wildcards.sample}.R1.duplicate_filter_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R2.fastq.gz".format(wildcards=wildcards)):
        out["R2"] = config.summary_directory + "{wildcards.sample}.R2.duplicate_filter_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.S.fastq.gz".format(wildcards=wildcards)):
        out["S"] = config.summary_directory + "{wildcards.sample}.S.duplicate_filter_summary.txt".format(wildcards=wildcards)
    return out

rule duplicate_filter_summary_combine1:
    input:
        unpack(duplicate_filter_summary_combine1_DetermineFiles)
    output:
        config.summary_directory + "{sample}.duplicate_filter_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        shell("head -1 %s > %s" %(input[0], output ) )
        for i in range(len(input)):
            shell("tail -n +2 %s >> %s" %(input[i], output))
        shell("rm -f {input}")

rule duplicate_filter_summary_combine2:
    input:
        expand( config.summary_directory + "{sample}.duplicate_filter_summary.txt", sample = SAMPLES )
    output:
        config.summary_directory + "duplicate_filter_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)
        shell("rm -f {input}")


rule quality_filter_paired:
    input:
        F=config.duplicate_filtered_directory + "{sample}.R1.dupfilter.fastq.gz",
        R=config.duplicate_filtered_directory + "{sample}.R2.dupfilter.fastq.gz"
    output:
        F=config.quality_filtered_directory + "{sample}.R1.fq.fastq.gz",
        R=config.quality_filtered_directory + "{sample}.R2.fq.fastq.gz",
        S=config.quality_filtered_directory + "{sample}.S.fq.temp2.fastq.gz"
    params:
        cluster="-cwd -l mfree=24G,h_rt=24:0:0"
    benchmark:
        "".join([config.log_directory, "quality_filter_paired.{sample}.log"])
    run:
        out_F_nonzip = output.F.rstrip(".gz")
        out_R_nonzip = output.R.rstrip(".gz")
        out_S_nonzip = output.S.rstrip(".gz")
        shell( " src/quality_filtering_wrapper.sh {input.F} --paired_fastq {input.R} {wildcards.sample} %s --paired_fastq_output %s --singleton_output %s" %(out_F_nonzip, out_R_nonzip, out_S_nonzip) )
        shell( "gzip %s" %(out_F_nonzip) )
        shell( "gzip %s" %(out_R_nonzip) )
        shell( "gzip %s" %(out_S_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.R}")
            shell("rm -f {input.F}")

rule quality_filter_singleton:
    input:
        S=config.duplicate_filtered_directory + "{sample}.S.dupfilter.fastq.gz"
    output:
        S=config.quality_filtered_directory + "{sample}.S.fq.temp.fastq.gz"
    params:
        cluster="-cwd -l mfree=24G,h_rt=24:0:0"
    benchmark:
        "".join([config.log_directory, "quality_filter_singleton.{sample}.log"])
    run:
        out_nonzip = output.S.rstrip(".gz")
        shell( "src/quality_filtering_wrapper.sh {input} {wildcards.sample} %s" %(out_nonzip) )
        shell( "gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.S}")

rule quality_filter_summary:
    input:
        pre_R1=config.duplicate_filtered_directory + "{sample}.R1.dupfilter.fastq.gz",
        pre_R2=config.duplicate_filtered_directory + "{sample}.R2.dupfilter.fastq.gz",
        post_R1=config.quality_filtered_directory + "{sample}.R1.fq.fastq.gz",
        post_R2=config.quality_filtered_directory + "{sample}.R2.fq.fastq.gz",
        new_singleton=config.quality_filtered_directory + "{sample}.S.fq.temp2.fastq.gz",
        filtered_singleton=config.quality_filtered_directory + "{sample}.S.fq.fastq.gz"
    output:
        config.summary_directory + "{sample}.quality_filter_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_quality_filtering.py {input.pre_R1} {input.pre_R2} {input.post_R1} {input.post_R2} {input.new_singleton} {input.filtered_singleton} > {output}"

rule quality_filter_summary_combine:
    input:
        expand( config.summary_directory + "{sample}.quality_filter_summary.txt", sample = SAMPLES )
    output:
        config.summary_directory + "quality_filter_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)


rule quality_filter_summary_script:
    input:
        config.quality_filtered_directory + "{sample}.{type}.fq.fastq.gz",
    output:
        config.summary_directory + "{sample}.{type}.fq.fastq_summary.txt",
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_input_fastq.py {input} > {output}"


def quality_filter_summary_basic_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.quality_filtered_directory + "{wildcards.sample}.R1.fq.fastq.gz".format(wildcards=wildcards)):
        out["R1"] = config.summary_directory + "{wildcards.sample}.R1.fq.fastq_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.quality_filtered_directory + "{wildcards.sample}.R2.fq.fastq.gz".format(wildcards=wildcards)):
        out["R2"] = config.summary_directory + "{wildcards.sample}.R2.fq.fastq_summary.txt".format(wildcards=wildcards)
    if os.path.isfile( config.quality_filtered_directory + "{wildcards.sample}.S.fq.fastq.gz".format(wildcards=wildcards)):
        out["S"] = config.summary_directory + "{wildcards.sample}.S.fq.fastq_summary.txt".format(wildcards=wildcards)
    return out

rule quality_filter_summary_basic:
    input:
        unpack(quality_filter_summary_basic_DetermineFiles)
    output:
        config.summary_directory + "{sample}.quality_filter_summary_basic.txt"
    params:
        cluster=default_cluster_params
    run:
        shell("head -1 %s > %s" %(input[0], output ) )
        for i in range(len(input)):
            shell("tail -n +2 %s >> %s" %(input[i], output))
        shell("rm -f {input}")


rule quality_filter_summary_basic_combine:
    input:
        expand( config.summary_directory + "{sample}.quality_filter_summary_basic.txt", sample = SAMPLES )
    output:
        config.summary_directory + "quality_filter_summary_basic.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)


def merge_singletons_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R1.fastq.gz".format(wildcards=wildcards)):
        out["R_S"] = config.quality_filtered_directory + "{wildcards.sample}.S.fq.temp2.fastq.gz".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.S.fastq.gz".format(wildcards=wildcards)):
        out["S"] = config.quality_filtered_directory + "{wildcards.sample}.S.fq.temp.fastq.gz".format(wildcards=wildcards)
    return out

rule merge_singletons:
    input:
        unpack(merge_singletons_DetermineFiles)
    output: 
        out=config.quality_filtered_directory + "{sample}.S.fq.fastq.gz"
    params:
        cluster=default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "zcat {input} > %s" %(out_nonzip) )
        shell( "gzip %s" %( out_nonzip ) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f %s" %(config.fastq_directory + "{wildcards.sample}.S.fq.temp2.fastq.gz") )
            shell("rm -f %s" %(config.fastq_directory + "{wildcards.sample}.S.fq.temp.fastq.gz") )


rule map_reads:
    input: config.quality_filtered_directory + "{sample}.{type}.fq.fastq.gz",
    output:
        zipped_output=config.diamond_output_directory + diamond_output_suffix + "/{sample}.{type}.gz"
    params:
        memory=config.memory,
        cpus=config.cpus,
        threads=config.cpus * 2,
        sensitivity=config.sensitivity,
        db=config.kegg_db_path + config.kegg_version + "/KEGG_" + config.kegg_version + "_" + config.taxon + ".dmnd",
        cluster = "-l mfree=%iG -l h_rt=24:00:00 -cwd -pe serial %i -q borenstein-short.q" %(config.memory / config.cpus,config.cpus)
    threads: config.cpus * 2
    benchmark:
        #"".join([config.log_directory, "map_reads.{sample}.{type}.", diamond_output_suffix, ".log"])
        config.log_directory + diamond_output_suffix + "/map_reads.{sample}.{type}.log"
    run:
        #Test if the fastq is empty
        c = 0
        print(input)
        with gzip.open(input[0]) as f:
            for line in f:
                c += 1
                if c > 1:
                    break
        if c == 0:
            shell( "touch %s" %(output.zipped_output.rstrip(".gz") ))
        else:
            shell( " ".join([ "/net/borenstein/vol1/PROGRAMS/diamond", "blastx", "--block-size", str(config.block_size), "--index-chunks", str(config.index_chunks), "--threads", str(params.threads), "--db", params.db, "--query", "{input}","--out", output.zipped_output.rstrip(".gz"), params.sensitivity ]) ),
        shell( " ".join([ "gzip", output.zipped_output.rstrip(".gz") ]) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

def combine_mapping_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.quality_filtered_directory + "{wildcards.sample}.R1.fq.fastq.gz".format(wildcards=wildcards)) | os.path.isfile( config.fastq_directory + "{wildcards.sample}.R1.fastq.gz".format(wildcards=wildcards)):
        out["R1"] = "".join([config.diamond_output_directory, diamond_output_suffix, "/{wildcards.sample}.R1.gz"]).format(wildcards=wildcards)
    if os.path.isfile( config.quality_filtered_directory + "{wildcards.sample}.R2.fq.fastq.gz".format(wildcards=wildcards)) | os.path.isfile( config.fastq_directory + "{wildcards.sample}.R2.fastq.gz".format(wildcards=wildcards)):
        out["R2"] = "".join([config.diamond_output_directory, diamond_output_suffix, "/{wildcards.sample}.R2.gz"]).format(wildcards=wildcards)
    if os.path.isfile( config.quality_filtered_directory + "{wildcards.sample}.S.fq.fastq.gz".format(wildcards=wildcards)) | os.path.isfile( config.fastq_directory + "{wildcards.sample}.S.fastq.gz".format(wildcards=wildcards)):
        out["S"] = "".join([config.diamond_output_directory, diamond_output_suffix, "/{wildcards.sample}.S.gz"]).format(wildcards=wildcards)
    return out

rule combine_mapping:
    input:
        unpack(combine_mapping_DetermineFiles)
    output:
        out=config.diamond_output_directory + diamond_output_suffix + "/{sample}.gz"
    params:
        cluster = default_cluster_params
    benchmark:
        #config.log_directory + diamond_output_suffix + "/{sample}.log"
        config.log_directory + "combine_mapping.{sample}.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "zcat {input} > %s" %( out_nonzip ) )
        shell( "gzip %s" %( out_nonzip ) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule map_reads_summary:
    input:
        config.diamond_output_directory + diamond_output_suffix + "/{sample}.gz"
    output:
        config.summary_directory + diamond_output_suffix + "/{sample}.map_reads_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_diamond.py {input} > {output}"

rule map_reads_summary_combine:
    input:
        expand( config.summary_directory + diamond_output_suffix + "/{sample}.map_reads_summary.txt", sample = SAMPLES )
    output:
        config.summary_directory + diamond_output_suffix + "/map_reads_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)


rule hit_filtering:
    input:
        config.diamond_output_directory + diamond_output_suffix + "/{sample}.gz"
    output:
        out=config.diamond_filtered_directory + hitfiltering_output_suffix + "/{sample}.diamond_filtered.gz"
    params:
        N=config.best_n_hits,
        filtering_method=config.filtering_method,
        kegg_version=config.filtering_method,
        cluster = default_cluster_params
    benchmark:
        "".join([config.log_directory, "diamond_filtered.{sample}.log"])
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell("src/filter_hits.py {input} {params.filtering_method} {params.kegg_version} -n {params.N} > %s " %(out_nonzip))
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule hit_filtering_summary:
    input:
        config.diamond_filtered_directory + hitfiltering_output_suffix + "/{sample}.diamond_filtered.gz"
    output:
        config.summary_directory + hitfiltering_output_suffix + "/{sample}.hit_filtering_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_hit_filtering.py {input} > {output}"

rule hit_filtering_summary_combine:
    input:
        expand( config.summary_directory + hitfiltering_output_suffix + "/{sample}.hit_filtering_summary.txt", sample = SAMPLES )
    output:
        config.summary_directory + hitfiltering_output_suffix + "/hit_filtering_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)


rule gene_mapper:
    input:
        config.diamond_filtered_directory + hitfiltering_output_suffix + "/{sample}.diamond_filtered.gz"
    output:
        out=config.diamond_counts_directory + genemapper_output_suffix + "/{sample}.genecounts.gz"
    params:
        count_method=config.count_method_gene,
        kegg_version=config.kegg_version,
        cluster=default_cluster_params
    benchmark:
        config.log_directory + "gene_mapper.{sample}.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell("src/count_genes.py {input} {wildcards.sample} {params.count_method} {params.kegg_version} --normalization length > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule gene_mapper_summary:
    input:
        config.diamond_counts_directory + genemapper_output_suffix + "/{sample}.genecounts.gz"
    output:
        config.summary_directory + genemapper_output_suffix + "/{sample}.gene_mapper_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_gene_counting.py {input} > {output}"

rule gene_mapper_summary_combine:
    input:
        expand( config.summary_directory + genemapper_output_suffix + "/{sample}.gene_mapper_summary.txt", sample = SAMPLES)
    output:
        config.summary_directory +  genemapper_output_suffix + "/gene_mapper_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)

rule ko_mapper:
    input:
        config.diamond_counts_directory + genemapper_output_suffix + "/{sample}.genecounts.gz"
    output:
        out=config.ko_counts_directory + komapper_output_suffix + "/{sample}.kocounts.gz"
    params:
        counting_method=config.count_method_ko,
        kegg_version=config.kegg_version,
        cluster=default_cluster_params
    benchmark:
        config.log_directory + "ko_mapper.{sample}.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/count_kos.py {input} {params.counting_method} {params.kegg_version} > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule ko_mapper_summary:
    input:
        config.ko_counts_directory + komapper_output_suffix + "/{sample}.kocounts.gz"
    output:
        config.summary_directory +  komapper_output_suffix + "/{sample}.ko_mapper_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_ko_counting.py {input} > {output}"

def combine_files(input, output):
    shell("head -1 %s > %s" %(input[0], output ) )
    for i in range(len(input)):
        shell("tail -n +2 %s >> %s" %(input[i], output))

rule ko_mapper_summary_combine:
    input:
        expand( config.summary_directory +  komapper_output_suffix + "/{sample}.ko_mapper_summary.txt", sample = SAMPLES)
    output:
        config.summary_directory +  komapper_output_suffix + "/ko_mapper_summary.txt"
    params:
        cluster=default_cluster_params
    run:
        combine_files(input, output)

rule merge_tables:
    input:
        expand( config.ko_counts_directory + komapper_output_suffix + "/{sample}.kocounts.gz", sample = SAMPLES)
    output:
        out=config.ko_counts_directory + komapper_output_suffix + "/merge_kocounts.gz"
    params:
        cluster=default_cluster_params
    benchmark:
        config.log_directory + "merge_tables.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/merge_tables.py {input} > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )

rule normalization:
    input:
        config.ko_counts_directory + komapper_output_suffix + "/merge_kocounts.gz"
    output:
        out=config.ko_normalized_directory + normalization_output_suffix + "/kocounts_normalized.gz"
    params:
        norm_method=config.norm_method,
        musicc_method=config.musicc_correction_method,
        cluster=default_cluster_params
    benchmark:
        config.log_directory + "normalization.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/normalization_wrapper.sh {input} {params.norm_method} %s --musicc_correct {params.musicc_method}"  %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )


rule ko_functional_summary:
    input:
        config.ko_normalized_directory + normalization_output_suffix + "/kocounts_normalized.gz"
    output:
        out=config.module_profiles_directory + functionalsummary_output_suffix + "/functionalsummary.gz"
    params:
        summary_method=config.summary_method,
        summary_level=config.summary_level,
        taxon=config.taxon,
        kegg_version=config.kegg_version,
        cluster=default_cluster_params,
    benchmark:
        config.log_directory + "ko_functional_summary.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/summarize_ko_to_higher_level_wrapper.sh {input} {params.summary_method} {params.summary_level} {params.taxon} {params.kegg_version} %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule functional_summary_summary: 
    input:
        config.module_profiles_directory + functionalsummary_output_suffix + "/functionalsummary.gz"
    output:
        #config.summary_directory + "functional_summary_summary." + functionalsummary_output_suffix + ".txt"
        config.summary_directory + functionalsummary_output_suffix + "/functional_summary_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_summarizing_ko_to_higher_level.py {input} > {output}"


def AllSummaries_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.fastq_directory + "%s.R1.fastq.gz" %(SAMPLES[0])) | os.path.isfile( config.fastq_directory + "%s.fq.fastq.gz" %(SAMPLES[0]) ) | os.path.isfile( config.fastq_directory + "%s.fq.fastq.gz" %(SAMPLES[0])):
        out["input_summary"] = config.summary_directory + "fastq_summary.txt"
        out["host_filtering_summary"] = config.summary_directory + "hostfilter_summary.txt"
        out["duplicate_filtering_summary"] = config.summary_directory + "duplicate_filter_summary.txt"
        out["quality_filtering_summary"] = config.summary_directory + "quality_filter_summary.txt"
    else:
        out["input_summary"] = config.summary_directory + "quality_filter_summary_basic.txt"
    out["mapping_summary"] = config.summary_directory + diamond_output_suffix + "/map_reads_summary.txt"
    out["blast_hit_filtering_summary"] = config.summary_directory + hitfiltering_output_suffix + "/hit_filtering_summary.txt"
    out["gene_counting_summary"] = config.summary_directory + genemapper_output_suffix + "/gene_mapper_summary.txt"
    out["ko_counting_summary"] = config.summary_directory +  komapper_output_suffix + "/ko_mapper_summary.txt"
    #out["functional_level_summarization_summaries"] = config.summary_directory + "functional_summary_summary." + functionalsummary_output_suffix + ".txt"
    out["functional_level_summarization_summaries"] = config.summary_directory + functionalsummary_output_suffix + "/functional_summary_summary.txt"
    return out

rule AllSummaries:
    input:
        unpack(AllSummaries_DetermineFiles)
    output:
        config.summary_directory + functionalsummary_output_suffix + "/Summaries_All.txt"
    params:
        cluster=default_cluster_params
    run:
        if os.path.isfile( config.fastq_directory + "%s.R1.fastq.gz" %(SAMPLES[0])) | os.path.isfile( config.fastq_directory + "%s.fq.fastq.gz" %(SAMPLES[0]) ) | os.path.isfile( config.fastq_directory + "%s.fq.fastq.gz" %(SAMPLES[0])): 
            shell("src/merge_pipeline_step_summary_tables.py {input.input_summary} --host_filtering_summary {input.host_filtering_summary} --duplicate_filtering_summary {input.duplicate_filtering_summary} --quality_filtering_summary {input.quality_filtering_summary} {input.mapping_summary} {input.blast_hit_filtering_summary} {input.gene_counting_summary} {input.ko_counting_summary} {input.functional_level_summarization_summaries} > {output}")
        else:
            shell("src/merge_pipeline_step_summary_tables.py {input.input_summary} {input.mapping_summary} {input.blast_hit_filtering_summary} {input.gene_counting_summary} {input.ko_counting_summary} {input.functional_level_summarization_summaries} > {output}")
