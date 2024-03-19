# Metagenomics workflow

Taxonomic classification of single reads from both amplicon-targeted and shotgun metagenomics sequencing.



## Introduction

This workflow can be used for the following:

+ Taxonomic classification of 16S rDNA and 18S rDNA amplicons using [default or custom databases](#FAQs). Default databases:
    - NCBI targeted loci: 16S rDNA, 18S rDNA, ITS (ncbi_16s_18s, ncbi_16s_18s_28s_ITS; see [here](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/) for details).
    - General databases: Standard-8, PlusPF-8, PlusPFP-8 (see [here](https://benlangmead.github.io/aws-indexes/k2) for details).
+ Generate taxonomic profiles of one or more metagenomic samples.
+ Identify [AMR genes](#4-identify-antimicrobial-resistance-genes-amr-optional).

Additional features:
+ Two different approaches are available: `kraken2` (k-mer based) or `minimap2` (using alignment).
+ Option to run it in [real time](#311-running-wf-metagenomics-in-real-time): `real_time`.
+ Results include:
    - An abundance table with counts per taxa in all the samples.
    - Interactive sankey and sunburst plots to explore the different identified lineages.
    - A bar plot comparing the abundances of the most abundant taxa in all the samples.




## Compute requirements

Recommended requirements:

+ CPUs = 12
+ Memory = 32GB

Minimum requirements:

+ CPUs = 6
+ Memory = 16GB

Approximate run time: ~40min for 1 million reads in total (24 barcodes) using Kraken2 and the Standard-8 database (using a previously downloaded db).

ARM processor support: True




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).  

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources. Therefore, nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of 
the required software. Both methods are automated out-of-the-box provided 
either Docker or Singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified in the example below. 

It is not required to clone or download the git repository in order to run the workflow. 
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository into the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-metagenomics --help 
```

A demo dataset is provided for testing of the workflow. It can be downloaded using: 

```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-metagenomics/wf-metagenomics-demo.tar.gz
tar -xzvf wf-metagenomics-demo.tar.gz
```

The workflow can be run with the demo data using: 

```
nextflow run epi2me-labs/wf-metagenomics \
--fastq wf-metagenomics-demo/test_data/ \
-profile standard 
```

For further information about running a workflow on the command line see https://labs.epi2me.io/wfquickstart/



## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| classifier | string | Kraken2 or Minimap2 workflow to be used for classification of reads. | Use Kraken2 for fast classification and minimap2 for finer resolution, see Readme for further info. | kraken2 |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| exclude_host | string | A FASTA or MMI file of the host reference. Reads that align with this reference will be excluded from the analysis. |  |  |


### Real Time Analysis Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| real_time | boolean | Enable to continuously watch the input directory for new input files. Reads will be analysed as they appear | This option enables the use of Nextflow’s directory watching feature to constantly monitor input directories for new files. As soon as files are written by an external process Nextflow will begin analysing these files. The workflow will accumulate data over time to produce an updating report. | False |
| batch_size | integer | Maximum number of sequence records to process in a batch. | Large files will be split such that batch_size records are processed together. Set to 0 to avoid rebatching input files. A value of 32000 is recommended to rebatch large files. | 0 |
| read_limit | integer | Stop processing data when a particular number of reads have been analysed. By default the workflow will run indefinitely. | Sets the upper bound on the number of reads that will be analysed before the workflow is automatically stopped and no more data is analysed. |  |
| port | integer | Network port for communication between Kraken2 server and clients (available in real time  pipeline). | The workflow uses a server process to handle Kraken2 classification requests. This allows the workflow to persist the sequence database in memory throughout the duration of processing. The option specifies the local network port on which the server and clients will communicate. | 8080 |
| host | string | Network hostname (or IP address) for communication between Kraken2 server and clients. (See also 'external_kraken2' parameter). (Available in real time  pipeline). | The workflow uses a server process to handle Kraken2 classification requests. This allows the workflow to persist the sequence database in memory throughout the duration of processing. The option specifies the local network hostname (or IP address) of the Kraken server. | localhost |
| external_kraken2 | boolean | Whether a pre-existing Kraken2 server should be used, rather than creating one as part of the workflow. (Available in real time  pipeline). | By default the workflow assumes that it is running on a single host computer, and further that it should start its own Kraken2 server. It may be desirable to start a Kraken2 server outside of the workflow, in which case this option should be enabled. This option may be used in conjunction with the `host` option to specify that the Kraken2 server is running on a remote computer.  | False |
| server_threads | integer | Number of CPU threads used by the Kraken2 server for classifying reads. (Available in the real_time pipeline). | For the real-time Kraken2 workflow, this is the number of CPU threads used by the Kraken2 server for classifying reads. | 2 |
| kraken_clients | integer | Number of clients that can connect at once to the Kraken-server for classifying reads. (Available in the real_time pipeline). | For the real-time Kraken2 workflow, this is the number of clients sending reads to the server. It should not be set to more than 4 fewer than the executor CPU limit. | 2 |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. Disabled in the real time pipeline. | The sample sheet is a CSV file with, minimally, columns named `barcode`,`alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. Disabled in the real time pipeline. |  |  |


### Reference Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| database_set | string | Sets the reference, databases and taxonomy datasets that will be used for classifying reads. Choices: ['ncbi_16s_18s','ncbi_16s_18s_28s_ITS', 'SILVA_138_1', 'Standard-8', 'PlusPF-8', 'PlusPFP-8']. Memory requirement will be slightly higher than the size of the database. Standard-8, PlusPF-8 and PlusPFP-8 databases require more than 8GB. | This setting is overridable by providing an explicit taxonomy, database or reference path in the other reference options. | Standard-8 |
| database | string | Not required but can be used to specifically override Kraken2 database [.tar.gz or Directory]. | By default uses database chosen in database_set parameter. |  |
| taxonomy | string | Not required but can be used to specifically override taxonomy database. Change the default to use a different taxonomy file  [.tar.gz or directory]. | By default NCBI taxonomy file will be downloaded and used. |  |
| reference | string | Override the FASTA reference file selected by the database_set parameter. It can be a FASTA format reference sequence collection or a minimap2 MMI format index. | This option should be used in conjunction with the database parameter to specify a custom database. |  |
| ref2taxid | string | Not required but can be used to specify a  ref2taxid mapping. Format is .tsv (refname  taxid), no header row. | By default uses ref2taxid for option chosen in database_set parameter. |  |
| taxonomic_rank | string | Returns results at the taxonomic rank chosen. In the Kraken2 pipeline: set the level that Bracken will estimate abundance at. Default: S (species). Other possible options are K (kingdom level), P (phylum), C (class), O (order), F (family), and G (genus). |  | S |


### Kraken2 Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bracken_length | integer | Set the length value Bracken will use | Should be set to the length used to generate the kmer distribution file supplied in the Kraken database input directory. For the default datasets these will be set automatically. ncbi_16s_18s = 1000 , ncbi_16s_18s_28s_ITS = 1000 , PlusPF-8 = 300 |  |
| kraken2_memory_mapping | boolean | Avoids loading database into RAM | Kraken 2 will by default load the database into process-local RAM; this flag will avoid doing so. It may be useful if the available RAM memory is lower than the size of the chosen database. | False |
| include_kraken2_assignments | boolean | A per sample TSV file that indicates how each input sequence was classified as well as the taxon that has been assigned to each read. The TSV's will only be output on completion of the workflow and therefore not at all if using the real time option whilst running indefinitely. |  | False |
| kraken2_confidence | number | Kraken2 Confidence score threshold. Default: 0.0. Valid interval: 0-1 | Apply a threshold to determine if a sequence is classified or unclassified. Please visit the following link for further details about how it works: https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring. | 0.0 |


### Minimap2 Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| minimap2filter | string | Filter output of minimap2 by taxids inc. child nodes, E.g. "9606,1404" | Provide a list of taxids if you are only interested in certain ones in your minimap2 analysis outputs. |  |
| minimap2exclude | boolean | Invert minimap2filter and exclude the given taxids instead | Exclude a list of taxids from analysis outputs. | False |
| split_prefix | boolean | Enable if using a very large reference with minimap2 | If reference fasta large enough to require multipart index, set to true to use split-prefix option with minimap2.  | False |
| keep_bam | boolean | Copy bam files into the output directory. |  | False |
| minimap2_by_reference | boolean | Add a table with the mean sequencing depth per reference, standard deviation and coefficient of variation. It adds a scatterplot of the sequencing depth vs. the coverage and a heatmap showing the depth per percentile to the report |  | False |


### Antimicrobial Resistance Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| amr | boolean | Scan reads for antimicrobial resistance or virulence genes | Reads will be scanned using abricate and the chosen database (`--amr_db`) to identify any acquired antimicrobial resistance or virulence genes found present in the dataset. NOTE: It cannot identify mutational resistance genes | False |
| amr_db | string | Database of antimicrobial resistance or virulence genes to use. |  | resfinder |
| amr_minid | integer | Threshold of required identity to report a match between a gene in the database and fastq reads. Valid interval: 0-100 |  | 80 |
| amr_mincov | integer | Minimum coverage (breadth-of) threshold required to report a match between a gene in the  database and fastq reads. Valid interval: 0-100 |  | 80 |


### Report Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| abundance_threshold | number | Remove those taxa whose abundance is equal or lower than the chosen value. | To remove taxa with abundances lower than or equal to a relative value (compared to the total number of reads), use a decimal between 0-1 (1 not inclusive). To remove taxa with abundances lower than or equal to an absolute value, provide a number larger than 1. | 0 |
| n_taxa_barplot | integer | Number of most abundant taxa to be displayed in the barplot. The rest of taxa will be grouped under the "Other" category. |  | 9 |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_len | integer | Specify read length lower limit. | Any reads shorter than this limit will not be included in the analysis. | 0 |
| min_read_qual | number | Specify read quality lower limit. | Any reads with a quality lower than this limit will not be included in the analysis. |  |
| max_len | integer | Specify read length upper limit | Any reads longer than this limit will not be included in the analysis. |  |
| threads | integer | Maximum number of CPU threads to use per workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the number of CPU threads for all such processes. The total CPU resource used by the workflow is constrained by the executor configuration. See server threads parameter for Kraken specific threads in the real_time pipeline. | 4 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-metagenomics-report.html | Report for all samples. | aggregated |
| Abundance table with counts per taxa | ./abundance_table_{{ taxonomic_rank }}.tsv | Per-taxa counts TSV, including all samples. | aggregated |
| Bracken report file | ./bracken/{{ alias }}.kraken2_bracken.report | TSV file with the abundance of each taxa. See more info here: https://github.com/jenniferlu717/Bracken#output-kraken-style-bracken-report. | per-sample |
| Kraken2 taxonomic assignment per read (Kraken2 pipeline) | ./kraken2/{{ alias }}.kraken2.report.txt | Lineage-aggregated counts. See more info here: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format. | per-sample |
| Kraken2 taxonomic asignment per read (Kraken2 pipeline) | ./kraken2/{{ alias }}.kraken2.assignments.tsv | TSV file with the taxonomic assignment per read. See more info here: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format. | per-sample |
| Host BAM file | ./host_bam/{{ alias }}.bam | BAM file generated from mapping filtered input reads to the host reference. | per-sample |
| BAM index file of host reads | ./host_bam/{{ alias }}.bai | BAM index file generated from mapping filtered input reads to the host reference. | per-sample |
| BAM file (minimap2) | ./bam/{{ alias }}.reference.bam | BAM file generated from mapping filtered input reads to the reference. | per-sample |
| BAM index file (minimap2) | ./bam/{{ alias }}.bai | Index file generated from mapping filtered input reads to the reference. | per-sample |
| JSON file with identified AMR genes (amr) | ./amr/{{ alias }}.json | JSON file with abricate results. See more info here: https://github.com/tseemann/abricate#output. | per-sample |




## Pipeline overview

### 1. Concatenate input files and generate per read stats

[fastcat](https://github.com/epi2me-labs/fastcat) is used to concatenate input FASTQ files prior to downstream processing of the workflow. It will also output per-read stats including read lengths and average qualities.

You may want to choose which reads are analysed by filtering them using these flags `max_len`, `min_len`, `min_read_qual`, (see the [Inputs section](#advanced-options) for details).

### 2. Remove host sequences (optional)

We have included an optional filtering step to remove any host sequences that map (using [Minimap2](https://github.com/lh3/minimap2)) against a provided host reference (e.g. human), which can be a FASTA file or a MMI index. To use this option provide the path to your host reference with the `exclude_host` parameter. The mapped reads are output in a BAM file and excluded from further analysis.

```
nextflow run epi2me-labs/wf-metagenomics --fastq test_data/case04/reads.fastq.gz --exclude_host test_data/case04/host.fasta.gz
```

### 3. Classify reads taxonomically

There are two different approaches to taxonomic classification:

#### 3.1 Using Kraken2

[Kraken2](https://github.com/DerrickWood/kraken2) provides the fastest method for the taxonomic classification of the reads. Then, [Bracken](https://github.com/jenniferlu717/Bracken) is used to provide an estimate of the species (or the selected taxonomic rank) abundance in the sample.

##### 3.1.1 Running wf-metagenomics in real time

The Kraken2 mode can be used in real-time, allowing the workflow to run parallel with an ongoing sequencing run as read data is being produced by the Oxford Nanopore Technologies sequencing instrument. In this case, [Kraken2](https://github.com/DerrickWood/kraken2) is used with the [Kraken2-server](https://github.com/epi2me-labs/kraken2-server) and the user can visualise the classification of reads and species abundances in a real-time updating report.    
In real-time mode, the workflow processes new input files as they become available in batches of the specified size. Thus, this option cannot be used with a single fastq as input.    
>Note: When using the workflow in real-time, the workflow will run indefinitely until a user interrupts the program (e.g with `ctrl+c` when on the command line). The workflow can be configured to complete automatically after a set number of reads have been analysed using the `read_limit` variable. Once this threshold has been reached, the program will emit a `STOP.fastq.gz` file into the fastq directory, which will instruct the workflow to complete. The "STOP.fastq.gz" file is then deleted.

```
nextflow run epi2me-labs/wf-metagenomics --fastq test_data/case01 --real_time --batch_size 1000 --read_limit 4000
```

If running the Kraken2 pipeline **real_time** in a cluster, there are two options to enable the workflow to be able to communicate with the Kraken-server: 

1. Run a Kraken-server separately outside of the workflow.
2. Submit the workflow job to run on a single node (i.e. running as if on a single local machine).

>Notes on CPU resource of Kraken-server and client in the real time workflow
The real-time subworkflow uses a server process to handle Kraken2 classification requests. This allows the workflow to persist the sequence database in memory throughout the duration of processing. There are some parameters that may be worth considering to improve the performance of the workflow:
+ port: The option specifies the local network port on which the server and clients will communicate.
+ host: Network hostname (or IP address) for communication between Kraken2 server and clients. (See also external_kraken2 parameter).
+ external_kraken2: Whether a pre-existing Kraken2 server should be used, rather than creating one as part of the workflow. By default the workflow assumes that it is running on a single host computer, and further that it should start its own Kraken2 server. It may be desirable to start a Kraken2 server outside of the workflow (for example to host a large database), in which case this option should be enabled. This option may be used in conjuction with the host option to specify that the Kraken2 server is running on a remote computer.
+ server_threads: Number of CPU threads used by the Kraken2 server for classifying reads.
+ kraken_clients: Number of clients that can connect at once to the Kraken-server for classifying reads. It should not be set to more than 4 fewer than the executor CPU limit.

#### 3.2 Using Minimap2

[Minimap2](https://github.com/lh3/minimap2) provides better resolution, but, depending on the reference database used, can take significantly more time. Also, running the workflow with minimap2 does not support real-time analysis.

```
nextflow run epi2me-labs/wf-metagenomics --fastq test_data/case01 --classifier minimap2
```

The creation of alignment statistics plots can be enabled with the `minimap2_by_reference` flag. Using this option produces a table and scatter plot in the report showing sequencing depth and coverage of each reference. The report also contains a heatmap indicating the sequencing depth over relative genomic coordinates for the references with the highest coverage (references with a mean coverage of less than 1% of the one with the largest value are omitted).

### 4. Identify Antimicrobial Resistance Genes (AMR) (optional)

The workflow can be used to determine the presence of acquired antimicrobial resistance (AMR) or virulence genes within the dataset. It uses [ABRicate](https://github.com/tseemann/abricate) to scan reads against a database of AMR/virulence genes.

```
nextflow run epi2me-labs/wf-metagenomics --fastq path/to/fastq/ --database_set PlusPF-8 --amr
```

>Note: ABRicate can only report the presence of acquired AMR/virulence genes but cannot identify SNP-mediated antimicrobial resistance. 

### 5. Prepare output

The main output of the wf-metagenomics pipeline is the `wf-metagenomics-report.html` which can be found in the output directory. It contains a summary of read statistics, the taxonomic composition of the sample and some diversity metrics. The results shown in the report can also be customised with several options. For example, you can use `abundance_threshold` to remove all taxa less prevalent than the threshold from the abundance table. When setting this parameter to a natural number, taxa with fewer absolute counts are removed. You can also pass a decimal between 0.0-1.0 to drop taxa of lower relative abundance. Furthermore, `n_taxa_barplot` controls the number of taxa displayed in the bar plot and groups the rest under the category ‘Other’.

The workflow output also contains Kraken and bracken reports for each sample. Additionally, the ‘species-abundance.tsv’ is a table with the counts of the different taxa per sample. You can use the flag `include_kraken2_assignments` to include a per sample TSV file that indicates how each input sequence was classified as well as the taxon that has been assigned to each read. This TSV file will only be output on completion of the workflow and therefore not at all if using the real time option whilst running indefinitely. This option is available in the Kraken2 pipeline.


#### 5.1 Diversity indices

Species diversity refers to the taxonomic composition in a specific microbial community. There are some useful concepts to take into account:
* Richness: number of unique taxonomic groups present in the community,
* Taxonomic group abundance: number of individuals of a particular taxonomic group present in the community,
* Evenness: refers to the equitability of the different taxonomic groups in terms of their abundances.
    Two different communities can host the same number of different taxonomic groups (i.e. they have the same richness), but they can have different evenness. For instance, if there is one taxon whose abundance is much larger in one community compared to the other.

There are three types of biodiversity measures described over a special scale <sup>[1](https://doi.org/10.2307/1218190), [2](https://doi.org/10.1016/B978-0-12-384719-5.00036-8)</sup>: alpha-, beta-, and gamma-diversity.
* Alpha-diversity refers to the richness that occurs within a community given area within a region.
* Beta-diversity defined as variation in the identities of species among sites, provides a direct link between biodiversity at local scales (alpha diversity) and the broader regional species pool (gamma diversity).
* Gamma-diversity is the total observed richness within an entire region.

To provide a quick overview of the alpha-diversity of the microbial community, we provide some of the most common diversity metrics calculated for a specific taxonomic rank <sup>[3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4224527/)</sup>, which can be chosen by the user with the `taxonomic_rank` parameter ('D'=Domain,'P'=Phylum, 'C'=Class, 'O'=Order, 'F'=Family, 'G'=Genus, 'S'=Species). By default, the rank is 'S' (species-level). Some of the included alpha diversity metrics are:

* Shannon Diversity Index (H): Shannon entropy approaches zero if a community is almost entirely made up of a single taxon.

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

* Berger-Parker dominance index (BP): expresses the proportional importance of the most abundant type, i.e., the ratio of number of individuals of most abundant species to the total number of individuals of all the species in the sample.

```math
BP = n_i/N
```
   where n<sub>i</sub> refers to the counts of the most abundant taxon and N is the total of counts.     


* Fisher’s alpha: Fisher (see Fisher, 1943<sup>[4](https://doi.org/10.2307/1411)</sup>) noticed that only a few species tend to be abundant while most are represented by only a few individuals ('rare biosphere'). These differences in species abundance can be incorporated into species diversity measurements such as the Fisher’s alpha. This index is based upon the logarithmic distribution of number of individuals of different species. 

```math
S = \alpha * ln(1 + N/\alpha)
```
   where S is the total number of taxa, N is the total number of individuals in the sample. The value of Fisher's $`\alpha`$ is calculated by iteration.

These indices are calculated by default using the original abundance table (see McMurdie and Holmes<sup>[5](https://pubmed.ncbi.nlm.nih.gov/24699258/)</sup>, 2014 and Willis<sup>[6](https://www.frontiersin.org/articles/10.3389/fmicb.2019.02407/full)</sup>, 2019). If you want to calculate them from a rarefied abundance table (i.e. all the samples have been subsampled to contain the same number of counts per sample, which is the 95% of the minimum number of total counts), you can download the rarefied table from the report.

The report also includes the rarefaction curve per sample which displays the mean of species richness for a subsample of reads (sample size). Generally, this curve initially grows rapidly, as most abundant species are sequenced and they add new taxa in the community, then slightly flattens due to the fact that 'rare' species are more difficult of being sampled, and because of that is more difficult to report an increase in the number of observed species.

> Note: Within each rank, each named taxon is a unique unit. The counts are the number of reads assigned to that taxon. All `Unknown` sequences are considered as a unique taxon



## Troubleshooting

+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).
+ When using the Minimap2 pipeline with a custom database, you must make sure that the `ref2taxid` and reference files are coherent, as well as the taxonomy database.
+ If your device doesn't have the resources to use large Kraken2 databases (e.g. Standard-8, PlusPF-8 and PlusPFP-8), you can enable `kraken2_memory_mapping` to reduce the amount of memory required.




## FAQ's

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-metagenomics/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 

+ *Which database is used by default?* - By default, the workflow uses the Standard-8 in kraken2 pipelines and the NCBI 16S + 18S rRNA database in the minimap2 workflow. It will be downloaded the first time the workflow is run and re-used in subsequent runs.

+ *Are more databases available?* - Other metagenomic databases (listed below) can be selected with the `database_set` parameter, but the workflow can also be used with a custom database if required (see [here](https://labs.epi2me.io/how-to-meta-offline/) for details).
    * 16S, 18S, ITS
        * ncbi_16s_18s and ncbi_16s_18s_28s_ITS:  Archaeal, bacterial and fungal 16S/18S and ITS data. There are two databases available using the data from [NCBI]https://www.ncbi.nlm.nih.gov/refseq/targetedloci/)
        * SILVA_138_1: The [SILVA](https://www.arb-silva.de/) database (version 138) is also available. Note that SILVA uses its own set of taxids, which do not match the NCBI taxids. We provide the respective taxdump files, but if you prefer using the NCBI ones, you can create them from the SILVA files ([NCBI](https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/ncbi/)). As the SILVA database uses genus level, the last taxonomic rank at which the analysis is carried out is genus (`taxonomic_rank G`).
    * General databases
        * Standard-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).
        * PlusPF-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core, protozoa and fungi. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).
        * PlusPFP-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core, protozoa, fungi and plant. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).

+ *How can I use Kraken2 indexes?* - There are different databases available [here](https://benlangmead.github.io/aws-indexes/k2).

+ *How can I use custom databases?* - If you want to run the workflow using your own Kraken2 database, you'll need to provide the database and an associated taxonomy dump. For a custom Minimap2 reference database, you'll need to provide a reference FASTA (or MMI) and an associated ref2taxid file. For a guide on how to build and use custom databases, take a look at our [article on how to run wf-metagenomics offline](https://labs.epi2me.io/how-to-meta-offline/).

+ *How can I run the workflow with less memory?* -
    When running in Kraken mode, you can set the `kraken2_memory_mapping` parameter if the available memory is smaller than the size of the database.

+ *How can I run the workflow offline?* - To run wf-metagenomics offline you can use the workflow to download the databases from the internet and prepare them for offline re-use later. If you want to use one of the databases supported out of the box by the workflow, you can run the workflow with your desired database and any input (for example, the test data). The database will be downloaded and prepared in a directory on your computer. Once the database has been prepared, it will be used automatically the next time you run the workflow without needing to be downloaded again. You can find advice on picking a suitable database in our [article on selecting databases for wf-metagenomics](https://labs.epi2me.io/metagenomic-databases/).

+ *Which databases are available for AMR?* - By default, ABRicate is set to search for AMR genes present in the [Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) database. Users can choose from a number of databases using the `amr_db` parameter. 

    |```amr_db``` | Database |
    |---------------|----------|
    |```resfinder```| [Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)|
    |```ecoli_vf```| [E. coli virulence factors](https://github.com/phac-nml/ecoli_vf)|
    |```plasmidfinder```| [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/)|
    |```card```| [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/)|
    |```argannot```| [ARG-ANNOT](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)|
    |```vfdb```| [Virulence factor DB](http://www.mgc.ac.cn/VFs/)|
    |```ncbi```| [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047)|
    |```megares```| [MEGAres](https://www.meglab.org/megares/)|
    |```ecoh```| [E. coli AMR DB from SRST2](https://github.com/katholt/srst2/tree/master/data)|



## Related blog posts

+ [How to build and use databases to run wf-metagenomics and wf-16s offline](https://labs.epi2me.io/how-to-meta-offline/).
+ [Selecting the correct databases in the wf-metagenomics](https://labs.epi2me.io/metagenomic-databases/).

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.




