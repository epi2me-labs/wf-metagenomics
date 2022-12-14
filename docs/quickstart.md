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

