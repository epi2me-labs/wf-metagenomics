# wf-metagenomics

wf-metagenomics is a Nextflow workflow for identification of the origin of single reads from both amplicon-targeted and shotgun metagenomics sequencing. The workflow has two modes of operation, it can use either [kraken2](https://ccb.jhu.edu/software/kraken2/) or [minimap2](https://github.com/lh3/minimap2) to determine the origin of reads.

The kraken2 mode can be used in real-time, allowing the workflow to run continuously alongside an ongoing sequencing run as read data is being produced by the Oxford Nanopore Technologies' sequencing instrument. The user can visualise the classification of reads and species abundances in a real-time updating report.

## Introduction

wf-metagenomics offers two different approaches to assigning sequence reads to a species:

### Kraken2 - Default

[Kraken2](https://github.com/DerrickWood/kraken2) is used with the [Kraken2-server](https://github.com/epi2me-labs/kraken2-server) to offer the fastest method for classification of reads. [Bracken](https://github.com/jenniferlu717/Bracken) is then used to give a good estimate of species level abundance in the sample which can be visualised in the report. The Kraken2 workflow mode can be run in real time. See quickstart below for more details.

### Minimap2 

[Minimap2](https://github.com/lh3/minimap2) provides the finest resolution analysis but, depending on the reference database used, at the expense of significantly more compute time. Currently the minimap2 mode does not support real-time.

The wf-metagenomics workflow by default uses the NCBI 16S + 18S rRNA database that will be downloaded at the start of an analysis, there are expanded metagenomic database options available with the --source parameter but the workflow is not tied to this database and can also be used with custom databases as required.

## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://sylabs.io/singularity/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-metagenomics --help
```

to see the options for the workflow.

The main options are 

* `fastq`: A fastq file or directory containing fastq input files or directories of input files. 
* `kraken2`: When set to true will run the analysis with Kraken2 and Bracken
* `minimap2`: When set to true will run the analysis with minimap2
* `watch_path`: Used to run the workflow in real-time, will continue to watch until a "STOP.fastq" is found
* `read_limit`: Used in combination with watch_path the specify an end point

***Kraken2***

You can run the workflow with test_data available in the github repository.

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data```

You can also run the workflow in real-time, meaning the workflow will watch the input directory(s) and process inputs at they become available in the batch sizes specified.

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data --watch_path --batch_size 1``` 

**Important Note**

When using the real-time functionality of the workflow, the input directory must contain sub-directories which themselves contain sequencing reads in fastq files. The is in contrast to the standard workflow which will accept reads provided as a single file or fastq files directly under the provided input directory.

The below is therefore the only input layout supported by the real-time functionality (the names of the child directories are unrestricted):

eg.

```
─── input_directory
    ├── barcode01
    │   ├── reads0.fastq
    │   └── reads1.fastq
    ├── barcode02
    │   ├── reads0.fastq
    │   ├── reads1.fastq
    │   └── reads2.fastq
    └── barcode03
        └── reads0.fastq
```

***Minimap2***

Alternatively you can run using minimap2 instead. Currently this mode does not support real-time.

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data --classifier minimap2```

**Workflow outputs**

The primary outputs of the workflow include:

* classified and unclassified reads,
* text files detailing lineages found,
* an HTML report document detailing the primary findings of the workflow.

## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)