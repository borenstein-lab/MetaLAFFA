"""
Pipeline operation parameters
-----------------------------

This configuration submodule contains parameters related to user operation of the pipeline.
"""

sample_list = None
"""
The file to specifies which samples to annotate. If "None", then samples will be determined based on the contents of the input FASTQ directory.
"""

delete_intermediates = False
"""
Whether output files from intermediate pipeline steps should be saved. If not, then they are deleted after they are no longer required.
"""