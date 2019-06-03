"""
Annotation parameters
---------------------

This configuration submodule contains parameters related to the steps of the annotation process.
"""

db_path = "/net/borenstein/vol1/DATA_DIAMONDDBs/"
"""
Path to the annotation databases. This likely won't need to be changed.
"""

db_type = "DIAMOND"
"""
Type of database to use. This may differ depending on the aligner used. For example, BLAST databases and DIAMOND databases are stored in different formats.

Currently supported options: ["DIAMOND"].
"""

db_source = "KEGG"
"""
Source of the annotated gene database (e.g., KEGG, MetaCyc). Available options may be aligner-specific depending on which formats each database has been generated in. 

Currently supported DIAMOND options: ["KEGG"].
"""

db_version = "2018_02_23"
"""
Version of the annotated gene database. This should normally indicate when the alignment database was generated.

Currently supported KEGG options: ["2013_07_15", "2018_02_23"].
"""

db_subtype = "prokaryotes"
"""
Subtype of the annotated gene database. Sometimes the database used (e.g., KEGG) should be generated to restrict searches to a specific subset of genes (e.g., only genes occurring in prokaryotic genomes). The exact subset required will depend on the purpose of the study. 

Currently supported KEGG 2013_08_15 options: ["bacteria"].

Currently supported KEGG 2018_02_23 options: ["prokaryotes"].
"""

aligner = "DIAMOND"
"""
Aligner. DIAMOND is currently the default because it is much faster than MBLASTX and the loss in sensitivity is negligible.

Currently supported options: ["DIAMOND"].
"""

alignment_options = "blastx --top 1 --evalue 0.001"
"""
Options for the aligner. These will likely need to be manually customized if using a non-default aligner because available options are aligner-specific. Defaults are set for DIAMOND with fastest speed, lowest sensitivity, and standard cutoffs for alignment hit quality.

Refer to the documentation page for an aligner to see available options.
"""

hit_filtering_method = "best_hit"
"""
Method for filtering out unwanted alignment hits. By default, we filter to only include hits with the best e-value for each read, which provides a high specificity, as demonstrated in {TODO citation to Rogan's paper}. Other options include best_ko, best_N_hits, or best_N_kos, which can improve sensitivity at the cost of specificity.

Current hit filtering options: ["best_hit", "best_ko", "best_N_hits", "best_N_kos"].
"""

hit_filtering_options = ""
"""
Options for the hit filtering method. These will likely need to be manually customized if using a non-default hit filtering method because available options are method-specific. Defaults are set for best_hit, which has no options.

Refer to the documentation page for a hit filtering method to see available options.
"""

gene_counting_method = "fractional"
"""
Method for calculating gene counts from alignments. By default, we use a fractional counting method that distributes the count for a single read among all genes the read mapped to after filtering hits. Another method, whole counting, is used by other annotation pipelines such as {TODO}, but this will double-count reads which may effect results by {TODO}.

Current gene counting methods: ["fractional", "whole"]
"""

ortholog_counting_method = "fractional"
"""
Method for calculating ortholog counts from gene counts. Similar to gene counting, by default, we use a fractional counting method that distributes the counts for a single gene among all orthology groups associated with that gene. Another method, whole counting, is used by other annotation pipelines such as {TODO}, but this will double-count genes which may effect results by {TODO}.

Current ortholog counting methods: ["fractional", "whole"]
"""

ortholog_normalization_method = "musicc"
"""
Method for normalizing ortholog counts. By default, we use MUSiCC {TODO citation}, which transforms ortholog counts into average copy number per genome based on universal single copy genes. This allows valid comparisons of ortholog abundances between different metagenomes.

Note that ortholog normalization methods are database source specific (e.g., MUSiCC uses KO ortholog annotations for normalization) and certain methods may not be available for some database sources.

Current KEGG ortholog normalization methods: ["musicc", "none"]
"""

ortholog_normalization_options = "--correct use_generic"
"""
Options for the ortholog normalization method. These will likely need to be manually customized if using a non-default normalization method because available options are method-specific. Defaults are set for MUSiCC with {TODO}.

Refer to the documentation page for an ortholog normalization method to see available options. 
"""

functional_summary_levels = "module,pathway"
"""
Functional levels for ortholog count summaries. This should consist of a comma-delimited list of available functional levels. These functional levels will depend on the database used (e.g., KEGG KOs can be summarized to modules, pathways, and functional levels in the BRITE hierarchy).

Current KEGG functional levels: ["module", "pathway"]
"""

functional_summary_method = "fractional"
"""
Method for calculating functional summary level counts from ortholog counts. Similar to gene counting, by default, we use a fractional counting method that distributes the counts for a single ortholog among all higher functional categories associated with that gene. Another method, whole counting, is used by other annotation pipelines such as {TODO}, but this will double-count genes which may effect results by {TODO}.

Current functional summary level counting methods: ["fractional", "whole"]
"""