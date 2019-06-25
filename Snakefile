# snakemake -p -c "qsub {params.cluster}" -j 50 --latency-wait 60

import config.operation as op
from parameter_parsing import *
from config.library_functions import *

# Configure parameters for each step
step_info = parse_pipeline_steps()
step_params = create_step_params(step_info)

# Get the list of samples
samples = generate_sample_list(step_info, step_params)

# Create dummy inputs if necessary
create_dummy_inputs(samples, step_info, step_params)

# Define restricted wildcard values from configuration module
types = op.wildcard_restrictions["type"]
mappings = op.wildcard_restrictions["mapping"]

# Create any missing output directories
make_directories()

rule all:
    input:
        step_params["final_outputs"]

rule fastq_summary:
    input:
        unpack(step_params["fastq_summary"]["input"])
    output:
        step_params["fastq_summary"]["output"]
    params:
        cluster=step_params["fastq_summary"]["cluster_params"]
    benchmark:
        step_params["fastq_summary"]["benchmark"]
    run:
        run_step(step_params["fastq_summary"], input, output, wildcards)

rule fastq_summary_combine:
    input:
        expand(list(step_params["fastq_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["fastq_summary_combine"]["output"]
    params:
        cluster=step_params["fastq_summary_combine"]["cluster_params"]
    benchmark:
        step_params["fastq_summary_combine"]["benchmark"]
    run:
        run_step(step_params["fastq_summary_combine"], input, output, wildcards)

rule host_filter:
    input:
        unpack(step_params["host_filter"]["input"])
    output:
        step_params["host_filter"]["output"]
    params:
        cluster=step_params["host_filter"]["cluster_params"]
    benchmark:
        step_params["host_filter"]["benchmark"]
    run:
        run_step(step_params["host_filter"], input, output, wildcards)

rule host_filter_summary:
    input:
        unpack(step_params["host_filter_summary"]["input"])
    output:
        step_params["host_filter_summary"]["output"]
    params:
        cluster=step_params["host_filter_summary"]["cluster_params"]
    benchmark:
        step_params["host_filter_summary"]["benchmark"]
    run:
        run_step(step_params["host_filter_summary"], input, output, wildcards)

rule host_filter_summary_combine:
    input:
        expand(list(step_params["host_filter_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["host_filter_summary_combine"]["output"]
    params:
        cluster=step_params["host_filter_summary_combine"]["cluster_params"]
    benchmark:
        step_params["host_filter_summary_combine"]["benchmark"]
    run:
        run_step(step_params["host_filter_summary_combine"], input, output, wildcards)

rule duplicate_filter:
    input:
        unpack(step_params["duplicate_filter"]["input"])
    output:
        step_params["duplicate_filter"]["output"]
    params:
        cluster=step_params["duplicate_filter"]["cluster_params"]
    benchmark:
        step_params["duplicate_filter"]["benchmark"]
    run:
        run_step(step_params["duplicate_filter"], input, output, wildcards)

rule duplicate_filter_summary:
    input:
        unpack(step_params["duplicate_filter_summary"]["input"])
    output:
        step_params["duplicate_filter_summary"]["output"]
    params:
        cluster=step_params["duplicate_filter_summary"]["cluster_params"]
    benchmark:
        step_params["duplicate_filter_summary"]["benchmark"]
    run:
        run_step(step_params["duplicate_filter_summary"], input, output, wildcards)

rule duplicate_filter_summary_combine:
    input:
        expand(list(step_params["duplicate_filter_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["duplicate_filter_summary_combine"]["output"]
    params:
        cluster=step_params["duplicate_filter_summary_combine"]["cluster_params"]
    benchmark:
        step_params["duplicate_filter_summary_combine"]["benchmark"]
    run:
        run_step(step_params["duplicate_filter_summary_combine"], input, output, wildcards)

rule quality_filter:
    input:
        unpack(step_params["quality_filter"]["input"])
    output:
        step_params["quality_filter"]["output"]
    params:
        cluster=step_params["quality_filter"]["cluster_params"]
    benchmark:
        step_params["quality_filter"]["benchmark"]
    run:
        run_step(step_params["quality_filter"], input, output, wildcards)

rule quality_filter_summary:
    input:
        unpack(step_params["quality_filter_summary"]["input"])
    output:
        step_params["quality_filter_summary"]["output"]
    params:
        cluster=step_params["quality_filter_summary"]["cluster_params"]
    benchmark:
        step_params["quality_filter_summary"]["benchmark"]
    run:
        run_step(step_params["quality_filter_summary"], input, output, wildcards)

rule quality_filter_summary_combine:
    input:
        expand(list(step_params["quality_filter_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["quality_filter_summary_combine"]["output"]
    params:
        cluster=step_params["quality_filter_summary_combine"]["cluster_params"]
    benchmark:
        step_params["quality_filter_summary_combine"]["benchmark"]
    run:
        run_step(step_params["quality_filter_summary_combine"], input, output, wildcards)

rule quality_filter_fastq_summary:
    input:
        unpack(step_params["quality_filter_fastq_summary"]["input"])
    output:
        step_params["quality_filter_fastq_summary"]["output"]
    params:
        cluster=step_params["quality_filter_fastq_summary"]["cluster_params"]
    benchmark:
        step_params["quality_filter_fastq_summary"]["benchmark"]
    run:
        run_step(step_params["quality_filter_fastq_summary"], input, output, wildcards)

rule quality_filter_fastq_summary_combine:
    input:
        expand(list(step_params["quality_filter_fastq_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["quality_filter_fastq_summary_combine"]["output"]
    params:
        cluster=step_params["quality_filter_fastq_summary_combine"]["cluster_params"]
    benchmark:
        step_params["quality_filter_fastq_summary_combine"]["benchmark"]
    run:
        run_step(step_params["quality_filter_fastq_summary_combine"], input, output, wildcards)

rule map_reads_to_genes:
    input:
        unpack(step_params["map_reads_to_genes"]["input"])
    output:
        step_params["map_reads_to_genes"]["output"]
    params:
        cluster=step_params["map_reads_to_genes"]["cluster_params"]
    threads:
        step_params["map_reads_to_genes"]["threads"]
    benchmark:
        step_params["map_reads_to_genes"]["benchmark"]
    run:
        run_step(step_params["map_reads_to_genes"], input, output, wildcards)

rule map_reads_to_genes_summary:
    input:
        unpack(step_params["map_reads_to_genes_summary"]["input"])
    output:
        step_params["map_reads_to_genes_summary"]["output"]
    params:
        cluster=step_params["map_reads_to_genes_summary"]["cluster_params"]
    benchmark:
        step_params["map_reads_to_genes_summary"]["benchmark"]
    run:
        run_step(step_params["map_reads_to_genes_summary"], input, output, wildcards)

rule map_reads_to_genes_summary_combine:
    input:
        expand(list(step_params["map_reads_to_genes_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["map_reads_to_genes_summary_combine"]["output"]
    params:
        cluster=step_params["map_reads_to_genes_summary_combine"]["cluster_params"]
    benchmark:
        step_params["map_reads_to_genes_summary_combine"]["benchmark"]
    run:
        run_step(step_params["map_reads_to_genes_summary_combine"], input, output, wildcards)

rule hit_filter:
    input:
        unpack(step_params["hit_filter"]["input"])
    output:
        step_params["hit_filter"]["output"]
    params:
        cluster=step_params["hit_filter"]["cluster_params"]
    benchmark:
        step_params["hit_filter"]["benchmark"]
    run:
        run_step(step_params["hit_filter"], input, output, wildcards)

rule hit_filter_summary:
    input:
        unpack(step_params["hit_filter_summary"]["input"])
    output:
        step_params["hit_filter_summary"]["output"]
    params:
        cluster=step_params["hit_filter_summary"]["cluster_params"]
    benchmark:
        step_params["hit_filter_summary"]["benchmark"]
    run:
        run_step(step_params["hit_filter_summary"], input, output, wildcards)

rule hit_filter_summary_combine:
    input:
        expand(list(step_params["hit_filter_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["hit_filter_summary_combine"]["output"]
    params:
        cluster=step_params["hit_filter_summary_combine"]["cluster_params"]
    benchmark:
        step_params["hit_filter_summary_combine"]["benchmark"]
    run:
        run_step(step_params["hit_filter_summary_combine"], input, output, wildcards)

rule hit_filter_combine:
    input:
        unpack(step_params["hit_filter_combine"]["input"])
    output:
        step_params["hit_filter_combine"]["output"]
    params:
        cluster=step_params["hit_filter_combine"]["cluster_params"]
    benchmark:
        step_params["hit_filter_combine"]["benchmark"]
    run:
        run_step(step_params["hit_filter_combine"], input, output, wildcards)

rule gene_map:
    input:
        unpack(step_params["gene_map"]["input"])
    output:
        step_params["gene_map"]["output"]
    params:
        cluster=step_params["gene_map"]["cluster_params"]
    benchmark:
        step_params["gene_map"]["benchmark"]
    run:
        run_step(step_params["gene_map"], input, output, wildcards)

rule gene_map_summary:
    input:
        unpack(step_params["gene_map_summary"]["input"])
    output:
        step_params["gene_map_summary"]["output"]
    params:
        cluster=step_params["gene_map_summary"]["cluster_params"]
    benchmark:
        step_params["gene_map_summary"]["benchmark"]
    run:
        run_step(step_params["gene_map_summary"], input, output, wildcards)

rule gene_map_summary_combine:
    input:
        expand(list(step_params["gene_map_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["gene_map_summary_combine"]["output"]
    params:
        cluster=step_params["gene_map_summary_combine"]["cluster_params"]
    benchmark:
        step_params["gene_map_summary_combine"]["benchmark"]
    run:
        run_step(step_params["gene_map_summary_combine"], input, output, wildcards)

rule ortholog_map:
    input:
        unpack(step_params["ortholog_map"]["input"])
    output:
        step_params["ortholog_map"]["output"]
    params:
        cluster=step_params["ortholog_map"]["cluster_params"]
    benchmark:
        step_params["ortholog_map"]["benchmark"]
    run:
        run_step(step_params["ortholog_map"], input, output, wildcards)

rule ortholog_map_summary:
    input:
        unpack(step_params["ortholog_map_summary"]["input"])
    output:
        step_params["ortholog_map_summary"]["output"]
    params:
        cluster=step_params["ortholog_map_summary"]["cluster_params"]
    benchmark:
        step_params["ortholog_map_summary"]["benchmark"]
    run:
        run_step(step_params["ortholog_map_summary"], input, output, wildcards)

rule ortholog_map_summary_combine:
    input:
        expand(list(step_params["ortholog_map_summary_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["ortholog_map_summary_combine"]["output"]
    params:
        cluster=step_params["ortholog_map_summary_combine"]["cluster_params"]
    benchmark:
        step_params["ortholog_map_summary_combine"]["benchmark"]
    run:
        run_step(step_params["ortholog_map_summary_combine"], input, output, wildcards)

rule ortholog_map_combine:
    input:
        expand(list(step_params["ortholog_map_combine"]["input"](None).values()), sample=samples, type=types)
    output:
        step_params["ortholog_map_combine"]["output"]
    params:
        cluster=step_params["ortholog_map_combine"]["cluster_params"]
    benchmark:
        step_params["ortholog_map_combine"]["benchmark"]
    run:
        run_step(step_params["ortholog_map_combine"], input, output, wildcards)

rule ortholog_abundance_correction:
    input:
        unpack(step_params["ortholog_abundance_correction"]["input"])
    output:
        step_params["ortholog_abundance_correction"]["output"]
    params:
        cluster=step_params["ortholog_abundance_correction"]["cluster_params"]
    benchmark:
        step_params["ortholog_abundance_correction"]["benchmark"]
    run:
        run_step(step_params["ortholog_abundance_correction"], input, output, wildcards)

rule ortholog_aggregation:
    input:
        unpack(step_params["ortholog_aggregation"]["input"])
    output:
        step_params["ortholog_aggregation"]["output"]
    params:
        cluster=step_params["ortholog_aggregation"]["cluster_params"]
    benchmark:
        step_params["ortholog_aggregation"]["benchmark"]
    run:
        run_step(step_params["ortholog_aggregation"], input, output, wildcards)

rule ortholog_aggregation_summary:
    input:
        unpack(step_params["ortholog_aggregation_summary"]["input"])
    output:
        step_params["ortholog_aggregation_summary"]["output"]
    params:
        cluster=step_params["ortholog_aggregation_summary"]["cluster_params"]
    benchmark:
        step_params["ortholog_aggregation_summary"]["benchmark"]
    run:
        run_step(step_params["ortholog_aggregation_summary"], input, output, wildcards)

rule ortholog_aggregation_summary_combine:
    input:
        expand(list(step_params["ortholog_aggregation_summary_combine"]["input"](None).values()), mapping=mappings)
    output:
        step_params["ortholog_aggregation_summary_combine"]["output"]
    params:
        cluster=step_params["ortholog_aggregation_summary_combine"]["cluster_params"]
    benchmark:
        step_params["ortholog_aggregation_summary_combine"]["benchmark"]
    run:
        run_step(step_params["ortholog_aggregation_summary_combine"], input, output, wildcards)

rule summary_combine:
    input:
        unpack(step_params["summary_combine"]["input"])
    output:
        step_params["summary_combine"]["output"]
    params:
        cluster=step_params["summary_combine"]["cluster_params"]
    benchmark:
        step_params["summary_combine"]["benchmark"]
    run:
        run_step(step_params["summary_combine"], input, output, wildcards)
