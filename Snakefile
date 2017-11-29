#Snakemake pipeline for Borenstein lab metagenomics annotation pipeline

# snakemake -p -c "qsub" -j 50

import config,gzip

diamond_output_suffix = "_".join([x for x in ["diamond", config.alignment_method, config.sensitivity, "top_percentage", str(config.top_percentage), "max_e_value", str(config.max_e_value)] if x])
#default_cluster_params = "-cwd -q borenstein-short.q"
default_cluster_params = "-cwd"

wildcard_constraints:
    sample="[A-Za-z0-9_]+",
    type="[A-Za-z0-9_]+"

rule all:
    input:
        #expand( config.output_directory + "{sample}.gz", sample = samples_all )
        #expand( config.fastq_directory + "{sample}.S.fq.fastq.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"] ),
        #expand( "".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]), sample = ["C68N1ACXX_7", "C82C3ACXX_1"], type = ["R1","R2","S"] )
        expand(config.diamond_output_directory + "{sample}_genecounts.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"])
        #expand(config.diamond_output_directory + "{sample}_diamond_filtered.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"])


#rule host_filter:

#rule duplicate_filter:

rule quality_filter_duplicate:
    input:
        F=config.fastq_directory + "{sample}.R1.fastq.gz",
        R=config.fastq_directory + "{sample}.R2.fastq.gz"
    output:
        F=config.fastq_directory + "{sample}.R1.fq.fastq.gz",
        R=config.fastq_directory + "{sample}.R2.fq.fastq.gz",
        S=config.fastq_directory + "{sample}.S.fq.temp2.fastq.gz"
    params:
        cluster=default_cluster_params
    run:
        out_F_nonzip = output.F.rstrip(".gz")
        out_R_nonzip = output.R.rstrip(".gz")
        out_S_nonzip = output.S.rstrip(".gz")
        shell( " src/quality_filtering_wrapper.sh {input.F} --paired_fastq {input.R} {wildcards.sample} %s --paired_fastq_output %s --singleton_output %s" %(out_F_nonzip, out_R_nonzip, out_S_nonzip) )
        shell( "gzip %s" %(out_F_nonzip) )
        shell( "gzip %s" %(out_R_nonzip) )
        shell( "gzip %s" %(out_S_nonzip) )
        #TODO, delete intermediate

rule quality_filter_single:
    input:
        S=config.fastq_directory + "{sample}.S.fastq.gz"
    output:
        S=config.fastq_directory + "{sample}.S.fq.temp.fastq.gz"
    params:
        cluster=default_cluster_params
    run:
        out_nonzip = output.S.rstrip(".gz")
        shell( "src/quality_filtering_wrapper.sh {input} {wildcards.sample} %s" %(out_nonzip) )
        shell( "gzip %s" %(out_nonzip) )
        #TODO, delete intermediate

def merge_singletons_DetermineFiles(wildcards):
    out = {}
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.R1.fastq.gz".format(wildcards=wildcards)):
        out["R_S"] = config.fastq_directory + "{wildcards.sample}.S.fq.temp2.fastq.gz".format(wildcards=wildcards)
    if os.path.isfile( config.fastq_directory + "{wildcards.sample}.S.fastq.gz".format(wildcards=wildcards)):
        out["S"] = config.fastq_directory + "{wildcards.sample}.S.fq.temp.fastq.gz".format(wildcards=wildcards)
    return out

rule merge_singletons:
    input:
        unpack(merge_singletons_DetermineFiles)
    output: 
        out=config.fastq_directory + "{sample}.S.fq.fastq.gz"
    params:
        cluster=default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "zcat {input} > %s" %(out_nonzip) )
        shell( "gzip %s" %( out_nonzip ) )
        #TODO, delete intermediate


rule map_reads:
    input: config.fastq_directory + "{sample}.{type}.fq.fastq.gz"
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

rule combine_mapping:
    input:
        lambda wildcards: expand( "".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]), sample = wildcards.sample, type = ["R1","R2","S"] )
    output:
        out=config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz"
    params:
        cluster = default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "zcat {input} > %s" %( out_nonzip ) )
        shell( "gzip %s" %( out_nonzip ) )

rule hit_filtering:
    input:
        config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz"
    output:
        out=config.diamond_output_directory + "{sample}_diamond_filtered.gz"
    params:
        N=config.best_n_hits,
        filtering_method=config.filtering_method,
        cluster = default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell("python src/filter_hits.py {input} {params.filtering_method} -n {params.N} > %s " %(out_nonzip))
        shell("gzip %s" %(out_nonzip) )

rule gene_mapper:
    input:
        config.diamond_output_directory + "{sample}_diamond_filtered.gz"
    output:
        out=config.diamond_output_directory + "{sample}_genecounts.gz"
    params:
        count_method=config.count_method,
        cluster = default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell("src/count_genes.py {input} {wildcards.sample} {params.count_method} --normalization length > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )

rule ko_mapper:
    input:
        config.diamond_output_directory + "{sample}_genecounts.gz"
    output:
        out=config.diamond_output_directory + "{sample}_kocounts.gz"
    params:
        counting_method=config.count_method,
        cluster = default_cluster_params
    run:
        out_nonzip = output.out.rstrip(".gz")
        shell( "src/count_kos.py {input} {params.counting_method} > %s" %(out_nonzip) )
        shell("gzip %s" %(out_nonzip) )

