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
* tsv files containing for a specific rank, the counts per taxa in each sample.

***Diversity***

Species diversity refers to the taxonomic composition in a specific microbial community. There are three main concepts:

* Richness: number of unique taxonomic groups  present in the community,
* Taxonomic group abundance: number of individuals of a particular taxonomic group present in the community,
* Evennes: refers to the equitability of the different taxonomic groups in terms of their abundances.

Two different communities can host the same number of different taxonomic groups (i.e. they have the same richness), but they can have different eveness if, for instance, in one of the communities there is one taxon whose abundance is much larger than the other' abundances.

To provide a quick overview of the diversity of the microbial community, we provide some of the most common indices calculated by a specific taxonomic rank <sup>[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4224527/)</sup>. This rank can be chosen by the user providind the flag *--bracken_level* and the desired rank: 'D'=Domain,'P'=Phylum, 'C'=Class, 'O'=Order, 'F'=Family, 'G'=Genus, 'S'=Species. By default, the rank is 'S' (species level). Some of these indices are:

* Shannon-Wiener Diversity Index (H): Shannon entropy approaches zero when one of the taxa is much more abundant than the others.    
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


These indices are calculated by default using the original abundance table (see McMurdie and Holmes<sup>[2](https://pubmed.ncbi.nlm.nih.gov/24699258/)</sup>, 2014 and Willis<sup>[3](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)</sup>, 2019). If you want to calculate them from a rarefied abundance table (i.e. all the samples have been subsampled to contain the same number of counts per sample, which is the 95% of the minimum number of total counts), you can use the *wf-metagenomics-rarefied.tsv* within the output directory.

The report also includes the rarefaction curve per sample which displays the mean of species richness for a subsample of reads (sample size). Generally, this curve initially grows rapidly, as most abundant species are sequenced and they add new taxa in the community, then slightly flattens due to the fact that 'rare' species are more difficult of being sampled, and because of that is more difficult to report an increase in the number of observed species.

*Note: Within each rank, each named taxon is considered to be an unique unit. The counts are the number of reads assigned to that taxon. All 'Unknown' sequences are considered as an unique taxon.*


