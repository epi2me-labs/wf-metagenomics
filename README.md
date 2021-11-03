# wf-metagenomics

This workflow can be used for both amplicon-targeted and shotgun metagenomics.

wf-metagenomics offers two different approaches to assigning sequence reads to a species:

- Kraken2 offers the fastest functionality and in combination with the accompanying bracken software can be used for a more quantitative assessment of taxonomic representation within a sample.
- Minimap2 provides the finest resolution analysis but, depending on the reference database used, at the expense of significantly more compute time.

The wf-metagenomics workflow by default uses the NCBI 16S + 18S rRNA database that will be downloaded at the start of an analysis. The workflow is not tied to this database and can also be used with custom databases as required.

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

**Workflow outputs**

The primary outputs of the workflow include:

* classified and unclassified reads,
* text files detailing lineages found,
* an HTML report document detailing the primary findings of the workflow.


## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)