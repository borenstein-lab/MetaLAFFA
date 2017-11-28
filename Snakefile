#Snakemake pipeline for Borenstein lab metagenomics annotation pipeline


import config,gzip
#from file_handling import *

samples_all = []
read_files_map = {}
read_files_all = []
#sample_file = custom_read(config.sample_file_name)
sample_file = open(config.sample_file_name)
line = sample_file.readline()
for line in sample_file:
    if line[0] != "#":
        line_split = line.strip().split('\t')
        if line_split[0] not in read_files_map:
            read_files_map[ line_split[0] ] = []
        read_files_map[ line_split[0] ] = [x.rstrip(".gz") for x in line_split[5:7] if x != "N/A"]
        read_files_all += read_files_map[ line_split[0] ]
        samples_all.append(line_split[0])

diamond_output_suffix = "_".join([x for x in ["diamond", config.alignment_method, config.sensitivity, "top_percentage", str(config.top_percentage), "max_e_value", str(config.max_e_value)] if x])


wildcard_constraints:
    sample="[A-Za-z0-9_]+",
    file="[A-Za-z0-9_]+",
    type="[A-Za-z0-9_]+"

rule all:
    input:
        #expand( config.output_directory + "{sample}.gz", sample = samples_all )
        expand( config.fastq_directory + "{sample}.S.fq.fastq.gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"] ),
        #expand( zipped_output="".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]), sample = ["C68N1ACXX_7", "C82C3ACXX_1"], type = ["R1","R2","S"] )
        expand( config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"] )
        #expand( config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz", sample = ["C68N1ACXX_7", "C82C3ACXX_1"] ),

#Dummy rule
rule head:
    input:
        lambda wilcards: expand(config.fastq_directory + "{file}.gz", file = read_files_map[wildcards.sample] )
    output:
        config.output_directory + "{file}_{sample}_head.gz"
    shell:
        "zcat {input} | head > {output}"

#Dummy rule
rule combine:
    input:
        lambda wildcards: expand( config.output_directory + "{file}_{sample}_head.gz", file = read_files_map[wildcards.sample] )
    output:
        config.output_directory + "{sample}.gz"
    shell:
        "zcat {input} > {output}"


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
    run:
        shell( " src/quality_filtering_wrapper.sh {input.F} --paired_fastq {input.R} {wildcards.sample} {output.F} --paired_fastq_output {output.R} --singleton_output {output.S}"  )

rule quality_filter_single:
    input:
        S=config.fastq_directory + "{sample}.S.fastq.gz"
    output:
        S=config.fastq_directory + "{sample}.S.fq.temp.fastq.gz"
    run:
        shell( "src/quality_filtering_wrapper.sh {input} {wildcards.sample} {output}")

def merge_singletons_DetermineFiles():
    if os.path.isfile(config.fastq_directory + "{wildcard.sample}.S.fq.temp.fastq.gz"):
        return config.fastq_directory + "{sample}.S.fq.temp.fastq.gz"
    else:
        return config.fastq_directory + "{sample}.S.fq.temp2.fastq.gz"

rule merge_singletons:
    input:
        merge_singletons_DetermineFiles()
    output: 
        config.fastq_directory + "{sample}.S.fq.fastq.gz"
    run:
        infile_1 = config.fastq_directory + "{wildcards.sample}.S.fq.temp.fastq.gz"
        infile_2 = config.fastq_directory + "{wildcards.sample}.S.fq.temp2.fastq.gz"
        shell( "touch {output}" )
        shell( "zcat %s >> {output" %(infile_1) )
        shell( "zact %s >> {output}" %(infile_2) )

#def map_reads_DetermineFiles(wildcards):
#    out = {}
#    if os.path.isfile(config.fastq_directory + "{wildcards.sample}.R1.fq.fastq.gz"):
#        out["R1"] = config.fastq_directory + "{sample}.R1.fq.fastq.gz"
#    if os.path.isfile(config.fastq_directory + "{wildcards.sample}.R2.fq.fastq.gz"):
#        out["R2"] = config.fastq_directory + "{sample}.R2.fq.fastq.gz"
#    if os.path.isfile(config.fastq_directory + "{wildcards.sample}.S.fq.fastq.gz"):
#        out["S"] = config.fastq_directory + "{sample}.S.fq.fastq.gz"
#    return out

#rule map_reads:
#    input: unpack(map_reads_DetermineFiles)
#    output:
#        #zipped_output="".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]),
#        config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz"
#    params:
#        memory=config.memory,
#        cpus=config.cpus,
#        threads=config.cpus * 2,
#        cluster = "-l mfree=15G -l h_rt=24:00:00 -cwd -pe serial 24 -q borenstein-short.q"
#    threads: config.cpus * 2
#    benchmark:
#        "".join([config.log_directory, "{sample}.", diamond_output_suffix, ".log"])
#    run:
#        print({input})
#        print("R1" in {input})
#        print({input}.R1)
#        qwe
#        for f in ["R1","R2","S"]:
#            if f in {input}:
#                shell( " ".join([ "/net/borenstein/vol1/PROGRAMS/diamond", "blastx", "--block-size", str(config.block_size), "--index-chunks", str(config.index_chunks), "--threads", str(params.threads), "--db", config.db, "--query", "{input}","--out", "{output.zipped_output}".rstrip(".gz")]) ),
#                shell( " ".join([ "gzip", "{output.zipped_output}".rstrip(".gz") ]) )

rule map_reads:
    input: config.fastq_directory + "{sample}.{type}.fq.fastq.gz"
    output:
        zipped_output="".join( [config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]),
    params:
        memory=config.memory,
        cpus=config.cpus,
        threads=config.cpus * 2,
        cluster = "-l mfree=15G -l h_rt=24:00:00 -cwd -pe serial 24 -q borenstein-short.q"
    threads: config.cpus * 2
    benchmark:
        "".join([config.log_directory, "{sample}.{type}.", diamond_output_suffix, ".log"])
    run:
        shell( " ".join([ "/net/borenstein/vol1/PROGRAMS/diamond", "blastx", "--block-size", str(config.block_size), "--index-chunks", str(config.index_chunks), "--threads", str(params.threads), "--db", config.db, "--query", "{input}","--out", "{output.zipped_output}".rstrip(".gz")]) ),
        shell( " ".join([ "gzip", "{output.zipped_output}".rstrip(".gz") ]) )

rule combine_mapping:
    input:
        lambda wildcards: expand( "".join([config.diamond_output_directory, "{sample}.{type}.",  diamond_output_suffix, ".gz"]), sample = wildcards.sample, type = ["R1","R2","S"] )
    output:
        config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz"
    shell:
        "zcat {input} > {output}"

rule hit_filtering:
    input:
        config.diamond_output_directory + "{sample}" + diamond_output_suffix + ".gz"
    output:
        config.diamond_output_directory + "{sample}_diamond_filtered.gz"
#    params:
#        N=XXX
#        filtering_method=XXX
#    shell: #GZIP?
#        "python src/filter_hits.py {input} {params.filtering_method} -n {params.N} > {output}"
    
#rule gene_mapper:

#rule ko_mapper:

#rule normalize:

#rule pathway_calculation:



#rule summarize_mapping:
#	input:
#		expand("".join([config.diamond_output_directory, "/{sample}.", diamond_output_suffix, ".gz"]), sample=samples)
#	output:
#		"mapping_summary.tab"
#	shell:
#		" ".join([config.summarize_mapping_wrapper,config.summarize_mapping,  "-o {output}", "{input}"])
