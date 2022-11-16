
## Introduction

wf-metagenomics offers two different approaches to assigning sequence reads to a species:

### Kraken2 

[Kraken2](https://github.com/DerrickWood/kraken2) is used with the [Kraken2-server](https://github.com/epi2me-labs/kraken2-server) to offer the fastest method for classification of reads. [Bracken](https://github.com/jenniferlu717/Bracken) is then used to give a good estimate of species level abundance in the sample which can be visualised in the report. The Kraken2 workflow mode can be run in real time. See quickstart below for more details.

### Minimap2 

[Minimap2](https://github.com/lh3/minimap2) provides the finest resolution analysis but, depending on the reference database used, at the expense of significantly more compute time. Currently the minimap2 mode does not support real-time.

The wf-metagenomics workflow by default uses the NCBI 16S + 18S rRNA database that will be downloaded at the start of an analysis, there are expanded metagenomic database options available with the --source parameter but the workflow is not tied to this database and can also be used with custom databases as required.

