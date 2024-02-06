+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ When using the Minimap2 pipeline with a custom database, you must make sure that the `ref2taxid` and reference files are coherent, as well as the taxonomy database.
+ If your device doesn't have the resources to use large Kraken2 databases (e.g. Standard-8, PlusPF-8 and PlusPFP-8), you can enable `kraken2_memory_mapping` to reduce the amount of memory required.
+ At this moment, the workflow does not support empty input files. Please make sure all the input files contain at least a few reads before starting the workflow and consider removing those empty file from the analysis. There can be empty barcode directories with no FASTQ files within them, but if the FASTQ file is there, it should contain reads.