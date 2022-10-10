# wf-metagenomics

This workflow can be used for both amplicon-targeted and shotgun metagenomics.

## Introduction

wf-metagenomics offers two different approaches to assigning sequence reads to a species:

### Kraken2 

[Kraken2](https://github.com/DerrickWood/kraken2) is used with the [Kraken2-server](https://github.com/epi2me-labs/kraken2-server) to offer the fastest method for classification of reads. [Bracken](https://github.com/jenniferlu717/Bracken) is then used to give a good estimate of species level abundance in the sample. This part of the workflow can be run in real time. See quickstart below for more details.

### Minimap2 

[Minimap2](https://github.com/lh3/minimap2) provides the finest resolution analysis but, depending on the reference database used, at the expense of significantly more compute time. Currently this method does not support real-time.

The wf-metagenomics workflow by default uses the NCBI 16S + 18S rRNA database that will be downloaded at the start of an analysis, there are expanded metagenomic database options available with the --source parameter but the workflow is not tied to this database and can also be used with custom databases as required.

## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-metagenomics --help
```

to see the options for the workflow.

The main options are 

*`fastq`: A fastq file or directory containing fastq input files or directories of input files. 
*`kraken2`: When set to true will run the analysis with Kraken2 and Bracken
*`minimap2`: When set to true will run the analysis with minimap2
*`watch_path`: Used to run the workflow in real-time, will continue to watch until a "STOP.fastq" is found
*`read_limit`: Used in combination with watch_path the specify an end point

***Kraken2***

You can run the workflow with test_data available in the github repository.

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data --kraken2```

You can also run the workflow in real-time, meaning the workflow will watch the input directory(s) and process inputs at they become available in the batch sizes specified.

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data --kraken2 --watch_path --batch_size 1``` 

*** Minimap2 ***

Alternatively you can run using minimap2 instead. Currently this mode does not support real-time.

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data --minimap2```

**Workflow outputs**

The primary outputs of the workflow include:

* classified and unclassified reads,
* text files detailing lineages found,
* an HTML report document detailing the primary findings of the workflow.

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)