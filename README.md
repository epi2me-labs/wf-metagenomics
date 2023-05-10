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

```nextflow run epi2me-labs/wf-metagenomics --fastq test_data --watch_path --batch_size 1000```

**Important Note**

When using the real-time functionality of the workflow, the input directory must contain sequencing reads in fastq files or sub-directories which themselves contain sequencing reads in fastq files. This is in contrast to the standard workflow which can additionally accept reads provided as a single file directly.

The below is therefore the only input layout supported by the real-time functionality (the names of the child directories are unrestricted):

eg.

```
 ─── input_directory        ─── input_directory
    ├── reads0.fastq            ├── barcode01
    └── reads1.fastq            │   ├── reads0.fastq
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

***Databases***

The wf-metagenomics pipeline has 4 pre-defined databases.

To analyze  archaeal, bacterial and fungal 16S/18S and ITS data, there are two databases available that we have put together using the data from [NCBI](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/). They can be used in both kraken2 and minimap2 pipelines:
* ncbi_16s_18s
* ncbi_16s_18s_28s_ITS

To analyze metagenomics data (not just 16S/18S rRNA and ITS) with the kraken2 pipeline, there are different databases available [here](https://benlangmead.github.io/aws-indexes/k2). We have selected two of them:
* PlusPF-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core, protozoa and fungi. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).
* PlusPFP-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core, protozoa, fungi and plant. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).

If you want to run the workflow using your own database, you can use the parameters: database_set, taxonomy, database (kraken2) and reference (either a FASTA format reference or a minimap2 MMI format index) and ref2taxid (minimap2). Run `nextflow run main.nf --help` to find out more about them.

***Output***

The main output of the wf-metagenomics pipeline is the `wf-metagenomics-report.html` which can be found in the output directory. It contains a summary of read statistics, the taxonomic composition of the community and some diversity metrics.

***Diversity***

Species diversity refers to the taxonomic composition in a specific microbial community. There are three main concepts:

* Richness: number of unique taxonomic groups present in the community,
* Taxonomic group abundance: number of individuals of a particular taxonomic group present in the community,
* Evenness: refers to the equitability of the different taxonomic groups in terms of their abundances.

Two different communities can host the same number of different taxonomic groups (i.e. they have the same richness), but they can have different evenness. For instance, if there is one taxon whose abundance is much larger in one community compared to the other.

To provide a quick overview of the diversity of the microbial community, we provide some of the most common indices calculated by a specific taxonomic rank <sup>[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4224527/)</sup>. This rank can be chosen by the user providind the flag *--bracken_level* and the desired rank: 'D'=Domain,'P'=Phylum, 'C'=Class, 'O'=Order, 'F'=Family, 'G'=Genus, 'S'=Species. By default, the rank is 'S' (species level). Some of these indices are:

* Shannon Diversity Index (H): Shannon entropy approaches zero when one of the taxa is much more abundant than the others.    
```math
H = -\sum_{i=1}^{S}p_i*ln(p_i)
```

* Simpson's Diversity Index (D): The range is from 0 (low diversity) to 1 (high diversity).    

```math
D = \sum_{i=1}^{S}p_i^2
```

* Pielou Index (J): The values range from 0 (presence of a dominant species) and 1 (maximum evennes).    

```math
J = H/ln(S)
```


These indices are calculated by default using the original abundance table (see McMurdie and Holmes<sup>[2](https://pubmed.ncbi.nlm.nih.gov/24699258/)</sup>, 2014 and Willis<sup>[3](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)</sup>, 2019). If you want to calculate them from a rarefied abundance table (i.e. all the samples have been subsampled to contain the same number of counts per sample, which is the 95% of the minimum number of total counts), you can use download the rarefied table from the report.

The report also includes the rarefaction curve per sample which displays the mean of species richness for a subsample of reads (sample size). Generally, this curve initially grows rapidly, as most abundant species are sequenced and they add new taxa in the community, then slightly flattens due to the fact that 'rare' species are more difficult of being sampled, and because of that is more difficult to report an increase in the number of observed species.

*Note: Within each rank, each named taxon is considered to be a unique unit. The counts are the number of reads assigned to that taxon. All 'Unknown' sequences are considered as a unique taxon.*






## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)