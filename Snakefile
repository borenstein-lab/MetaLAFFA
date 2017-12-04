# Snakemake pipeline for Borenstein lab metagenomics annotation pipeline
# Author: Adrian Verster
# December 2017

# snakemake -p -c "qsub {params.cluster}" -j 50

import config,gzip,os

diamond_output_suffix = "_".join([x for x in ["diamond", config.alignment_method, config.sensitivity, "top_percentage", str(config.top_percentage), "max_e_value", str(config.max_e_value)] if x])
default_cluster_params = "-cwd -l mfree=20G"
delete_intermediates = False

#Without this line snakemake will sometimes fail a job because it fails to detect the output file due to latency
shell.suffix("; sleep 40")

#SAMPLES = ["C68N1ACXX_7", "C82C3ACXX_1"]
SAMPLES = list(set([x.split(".")[0] for x in os.listdir(config.fastq_directory)]))

wildcard_constraints:
    sample="[A-Za-z0-9_]+",
    type="[A-Za-z0-9_]+"

rule all:
    input:
        #expand( config.output_directory + "{sample}.gz", sample = samples_all )
        #expand( config.fastq_directory + "{sample}.S.fq.fastq.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"] ),
        #expand( "".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]), sample = ["C68N1ACXX_7", "C82C3ACXX_1"], type = ["R1","R2","S"] )
        #expand(config.diamond_output_directory + "{sample}_genecounts.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"]),
        #expand(config.diamond_output_directory + "{sample}_kocounts_normalized.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"])
        config.module_profiles_directory + "functionalsummary.gz",
        #Summaries
        expand(config.summary_directory + "{sample}.{type}.host_filer_summary.txt", sample = SAMPLES, type = ["R1","R2","S"]),
        expand(config.summary_directory + "{sample}.map_reads_summary.txt", sample = SAMPLES),
        expand(config.summary_directory + "{sample}.hit_filtering_summary.txt", sample = SAMPLES),
        expand(config.summary_directory + "{sample}.gene_mapper_summary.txt", sample = SAMPLES),
        expand(config.summary_directory + "functional_summary_summary.txt", sample = SAMPLES),



rule host_filter_duplicate:
    input:
        R1=config.fastq_directory + "{sample}.R1.fastq.gz",
        R2=config.fastq_directory + "{sample}.R2.fastq.gz"
    output:
        R1=config.host_filtered_directory + "{sample}.R1.hostfilter.fastq.gz",
        R2=config.host_filtered_directory + "{sample}.R2.hostfilter.fastq.gz",
        host_marked=config.host_filtered_directory + "{sample}.hostmarked_dup.fastq.gz"
    params:
        cluster=default_cluster_params
    run:
        out_F_nonzip = output.R1.rstrip(".gz")
        out_R_nonzip = output.R2.rstrip(".gz")
        out_host_nonzip = output.host_marked.rstrip(".gz")
        shell( "src/host_filtering_wrapper.sh {input.R1} %s %s --paired_fastq {input.R2} --paired_fastq_output %s" %(out_host_nonzip, out_F_nonzip, out_R_nonzip) )
        shell( "gzip %s" %(out_F_nonzip) )
        shell( "gzip %s" %(out_R_nonzip) )
        shell( "gzip %s" %(out_host_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.R1}")
            shell("rm -f {input.R2}")

rule host_filter_singleton:
    input:
        S=config.fastq_directory + "{sample}.S.fastq.gz",
    output:
        S=config.host_filtered_directory + "{sample}.S.hostfilter.fastq.gz",
        host_marked=config.host_filtered_directory + "{sample}.hostmarked_sig.fastq.gz"
    params:
        cluster=default_cluster_params
    run:
        out_nonzip = output.S.rstrip(".gz")
        out_host_nonzip = output.host_marked.rstrip(".gz")
        shell(" src/host_filtering_wrapper.sh {input.S} %s %s" %(out_host_nonzip, out_nonzip) )
        shell( "gzip %s" %(out_nonzip) )
        shell( "gzip %s" %(out_host_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.S}")

rule host_filer_summary:
    input:
        config.host_filtered_directory + "{sample}.{type}.hostfilter.fastq.gz"
    output:
        config.summary_directory + "{sample}.{type}.host_filer_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_host_filtering.py {input} > {output}"

rule duplicate_filter_singleton:
    input:
        S=config.host_filtered_directory + "{sample}.S.hostfilter.fastq.gz",
    output:
        S=config.duplicate_filtered_directory + "{sample}.S.dupfilter.fastq.gz",
        metrics_output=config.duplicate_filtered_directory + "{sample}.S.metricout.fastq.gz",
        duplicate_marked_reads=config.duplicate_filtered_directory + "{sample}.S.duplicatemarked.fastq.gz",
    params:
        cluster=default_cluster_params
    run:
        out_nonzip = output.S.rstrip(".gz")
        out_metrics_nonzip = output.metrics_output.rstrip(".gz")
        out_duplicate_nonzip = output.duplicate_marked.rstrip(".gz")
        shell(" src/duplicate_filtering_wrapper.sh {input.S} {wildcards.sample} %s %s %s " %(out_metrics_nonzip, out_duplicate_nonzip, out_nonzip) )
        shell( "gzip %s" %(out_S_nonzip) )
        shell( "gzip %s" %(out_metrics_nonzip) )
        shell( "gzip %s" %(out_duplicate_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.S}")


rule duplicate_filter_paired:
    input:
        R1=config.host_filtered_directory + "{sample}.R1.hostfilter.fastq.gz",
        R2=config.host_filtered_directory + "{sample}.R2.hostfilter.fastq.gz",
    output:
        R1=config.duplicate_filtered_directory + "{sample}.R1.dupfilter.fastq.gz",
        R2=config.duplicate_filtered_directory + "{sample}.R2.dupfilter.fastq.gz",
    params:
        metrics_output="",
        duplicate_marked_reads="",
        cluster=default_cluster_params
    run:
        out_F_nonzip = output.F.rstrip(".gz")
        out_R_nonzip = output.R.rstrip(".gz")
        shell( "src/duplicate_filtering_wrapper.sh {input.R1} {wildcard.sample} {params.metrics_output} {params.duplicate_marked_reads} %s --paired_fastq {input.R2} --paired_fastq_output %s" %(out_F_nonzip, out_R_nonzip) )
        shell( "gzip %s" %(out_F_nonzip) )
        shell( "gzip %s" %(out_R_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.R}")
            shell("rm -f {input.F}")

rule duplicate_filter_summary:
    input:
        config.duplicate_filtered_directory + "{sample}.{type}.dupfilter.fastq.gz"
    output:
        config.summary_directory + "{sample}.{type}.duplicate_filer_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "#src/summarize_duplicate_filtering.py {input} > {output}"

rule quality_filter_duplicate:
    input:
        F=config.duplicate_filtered_directory + "{sample}.R1.dupfilter.fastq.gz",
        R=config.duplicate_filtered_directory + "{sample}.R2.dupfilter.fastq.gz"
    output:
        F=config.quality_filtered_directory + "{sample}.R1.fq.fastq.gz",
        R=config.quality_filtered_directory + "{sample}.R2.fq.fastq.gz",
        S=config.quality_filtered_directory + "{sample}.S.fq.temp2.fastq.gz"
    params:
        cluster="-cwd -l mfree=24G"
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

rule quality_filter_single:
    input:
        S=config.duplicate_filtered_directory + "{sample}.S.dupfilter.fastq.gz"
    output:
        S=config.quality_filtered_directory + "{sample}.S.fq.temp.fastq.gz"
    params:
        cluster="-cwd -l mfree=24G"
    run:
        out_nonzip = output.S.rstrip(".gz")
        shell( "src/quality_filtering_wrapper.sh {input} {wildcards.sample} %s" %(out_nonzip) )
        shell( "gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input.S}")

rule quality_filter_summary:
    input:
        pre_R1=config.quality_filtered_directory + "{sample}.R1.dupfilter.fastq.gz",
        pre_R2=config.quality_filtered_directory + "{sample}.R2.dupfilter.fastq.gz",
        post_R1=config.fastq_directory + "{sample}.R1.fq.fastq.gz",
        post_R2=config.fastq_directory + "{sample}.R2.fq.fastq.gz",
        new_singleton=config.quality_filtered_directory + "{sample}.S.fq.temp2.fastq.gz",
        filtered_singleton=config.quality_filtered_directory + "{sample}.S.fq.fastq.gz"
    output:
        config.summary_directory + "{sample}.{type}.quality_filer_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_quality_filtering.py {input.pre_R1} {input.pre_R2} {input.post_R1} {input.post_R2} {input.new_singleton} {input.filtered_singleton} > {output}"



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
    input: config.quality_filtered_directory + "{sample}.{type}.fq.fastq.gz"
    output:
        zipped_output="".join( [config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]),
    params:
        memory=config.memory,
        cpus=config.cpus,
        threads=config.cpus * 2,
        cluster = "-l mfree=10G -l h_rt=24:00:00 -cwd -pe serial 24 -q borenstein-short.q"
    threads: config.cpus * 2
    benchmark:
        "".join([config.log_directory, "{sample}.{type}.", diamond_output_suffix, ".log"])
    run:
        shell( " ".join([ "/net/borenstein/vol1/PROGRAMS/diamond", "blastx", "--block-size", str(config.block_size), "--index-chunks", str(config.index_chunks), "--threads", str(params.threads), "--db", config.db, "--query", "{input}","--out", output.zipped_output.rstrip(".gz")]) ),
        shell( " ".join([ "gzip", output.zipped_output.rstrip(".gz") ]) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )


rule combine_mapping:
    input:
        lambda wildcards: expand( "".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]), sample = wildcards.sample, type = ["R1","R2","S"] )
    output:
        out=config.diamond_output_directory + "{sample}." + diamond_output_suffix + ".gz"
    params:
        cluster = default_cluster_params
    benchmark:
        config.log_directory + "{sample}" + diamond_output_suffix + ".log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "zcat {input} > %s" %( out_nonzip ) )
        shell( "gzip %s" %( out_nonzip ) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule map_reads_summary:
    input:
        config.diamond_output_directory + "{sample}." + diamond_output_suffix + ".gz"
    output:
        config.summary_directory + "{sample}.map_reads_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_diamond.py {input} > {output}"

rule hit_filtering:
    input:
        config.diamond_output_directory + "{sample}." + diamond_output_suffix + ".gz"
    output:
        out=config.diamond_filtered_directory + "{sample}_diamond_filtered.gz"
    params:
        N=config.best_n_hits,
        filtering_method=config.filtering_method,
        cluster = default_cluster_params
    benchmark:
        config.log_directory + "{sample}_diamond_filtered.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell("src/filter_hits.py {input} {params.filtering_method} -n {params.N} > %s " %(out_nonzip))
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule hit_filtering_summary:
    input:
        config.diamond_filtered_directory + "{sample}_diamond_filtered.gz"
    output:
        config.summary_directory + "{sample}.hit_filtering_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_hit_filtering.py {input} > {output}"

rule gene_mapper:
    input:
        config.diamond_filtered_directory + "{sample}_diamond_filtered.gz"
    output:
        out=config.gene_counts_directory + "{sample}_genecounts.gz"
    params:
        count_method=config.count_method,
        cluster = "-cwd -l mfree=2G"
    benchmark:
        config.log_directory + "{sample}_genecounts.log"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell("src/count_genes.py {input} {wildcards.sample} {params.count_method} --normalization length > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule gene_mapper_summary:
    input:
        config.gene_counts_directory + "{sample}_genecounts.gz"
    output:
        config.summary_directory + "{sample}.gene_mapper_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_gene_counting.py {input} > {output}"

rule ko_mapper:
    input:
        config.gene_counts_directory + "{sample}_genecounts.gz"
    output:
        out=config.ko_counts_directory + "{sample}_kocounts.gz"
    params:
        counting_method=config.count_method,
        cluster = "-cwd -l mfree=4G"
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/count_kos.py {input} {params.counting_method} > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule ko_mapper_summary:
    input:
        config.ko_counts_directory + "{sample}_kocounts.gz"
    output:
        config.summary_directory + "{sample}.ko_mapper_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_ko_counting.py {single_sample_ko_abundances} > {output}"

rule merge_tables:
    input:
        expand( config.ko_counts_directory + "{sample}_kocounts.gz", sample = SAMPLES)
    output:
        out=config.ko_counts_directory + "merge_kocounts.gz"
    params:
        cluster=default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/merge_tables.py {input} > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )

rule normalization:
    input:
        config.ko_counts_directory + "merge_kocounts.gz"
    output:
        out=config.ko_normalized_directory + "kocounts_normalized.gz"
    params:
        norm_method=config.norm_method,
        musicc_method=config.musicc_correction_method,
        cluster=default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/normalization_wrapper.sh {input} {params.norm_method} %s --musicc_correct {params.musicc_method}"  %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )


rule ko_functional_summary:
    input:
        config.ko_normalized_directory + "kocounts_normalized.gz"
    output:
        out=config.module_profiles_directory + "functionalsummary.gz"
    params:
        mapping_matrix=config.mapping_matrix,
        summary_method=config.summary_method,
        cluster=default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/summarize_ko_to_higher_level_wrapper.sh {input} {params.summary_method} {params.mapping_matrix} %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )
        #Delete intermediate
        if delete_intermediates:
            shell("rm -f {input}" )

rule functional_summary_summary: 
    input:
        config.module_profiles_directory + "functionalsummary.gz"
    output:
        config.summary_directory + "functional_summary_summary.txt"
    params:
        cluster=default_cluster_params
    shell:
        "src/summarize_summarizing_ko_to_higher_level.py {input} > {output}"
