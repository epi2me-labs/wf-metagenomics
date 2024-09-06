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

In addition, the user can output BAM files in a folder called `bams` by using the option `keep_bam`. If the user provides a custom database and uses the `igv` option, the workflow will also output the references with reads mappings, as well as an IGV configuration file. This configuration file allows the user to view the alignments in the EPI2ME Desktop Application in the Viewer tab. Note that the number of references can be reduced using the `abundance_threshold` option, which will select those references with a number of reads aligned higher than this value. Please, consider that the view of the alignment is highly dependent on the reference selected.

### 4. Identify Antimicrobial Resistance Genes (AMR) (optional)

The workflow can be used to determine the presence of acquired antimicrobial resistance (AMR) or virulence genes within the dataset. It uses [ABRicate](https://github.com/tseemann/abricate) to scan reads against a database of AMR/virulence genes.

```
nextflow run epi2me-labs/wf-metagenomics --fastq path/to/fastq/ --database_set PlusPF-8 --amr
```

>Note: ABRicate can only report the presence of acquired AMR/virulence genes but cannot identify SNP-mediated antimicrobial resistance. 

### 5. Prepare output

The main output of the wf-metagenomics pipeline is the `wf-metagenomics-report.html` which can be found in the output directory. It contains a summary of read statistics, the taxonomic composition of the sample and some diversity metrics. The results shown in the report can also be customised with several options. For example, you can use `abundance_threshold` to remove all taxa less prevalent than the threshold from the abundance table. When setting this parameter to a natural number, taxa with fewer absolute counts are removed. You can also pass a decimal between 0.0-1.0 to drop taxa of lower relative abundance. Furthermore, `n_taxa_barplot` controls the number of taxa displayed in the bar plot and groups the rest under the category ‘Other’.

The workflow output also contains Kraken and bracken reports for each sample. Additionally, the ‘species-abundance.tsv’ is a table with the counts of the different taxa per sample.


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