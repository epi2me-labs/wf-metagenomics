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
| database_set | string | Sets the reference, databases and taxonomy datasets that will be used for classifying reads. Choices: ['ncbi_16s_18s','ncbi_16s_18s_28s_ITS', 'SILVA_138_1', 'Standard-8', 'PlusPF-8', 'PlusPFP-8']. Memory requirement will be slightly higher than the size of the database. Standard-8, PlusPF-8 and PlusPFP-8 databases require more than 8GB and are only available in the kraken2 approach. | This setting is overridable by providing an explicit taxonomy, database or reference path in the other reference options. | Standard-8 |
| database | string | Not required but can be used to specifically override Kraken2 database [.tar.gz or Directory]. | By default uses database chosen in database_set parameter. |  |
| taxonomy | string | Not required but can be used to specifically override taxonomy database. Change the default to use a different taxonomy file  [.tar.gz or directory]. | By default NCBI taxonomy file will be downloaded and used. |  |
| reference | string | Override the FASTA reference file selected by the database_set parameter. It can be a FASTA format reference sequence collection or a minimap2 MMI format index. | This option should be used in conjunction with the database parameter to specify a custom database. |  |
| ref2taxid | string | Not required but can be used to specify a  ref2taxid mapping. Format is .tsv (refname  taxid), no header row. | By default uses ref2taxid for option chosen in database_set parameter. |  |
| taxonomic_rank | string | Returns results at the taxonomic rank chosen. In the Kraken2 pipeline: set the level that Bracken will estimate abundance at. Default: S (species). Other possible options are P (phylum), C (class), O (order), F (family), and G (genus). |  | S |


### Kraken2 Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| bracken_length | integer | Set the length value Bracken will use | Should be set to the length used to generate the kmer distribution file supplied in the Kraken database input directory. For the default datasets these will be set automatically. ncbi_16s_18s = 1000 , ncbi_16s_18s_28s_ITS = 1000 , PlusPF-8 = 300 |  |
| bracken_threshold | integer | Set the minimum read threshold Bracken will use to consider a taxon | Bracken will only consider taxa with a read count greater than or equal to this value. | 10 |
| kraken2_memory_mapping | boolean | Avoids loading database into RAM | Kraken 2 will by default load the database into process-local RAM; this flag will avoid doing so. It may be useful if the available RAM memory is lower than the size of the chosen database. | False |
| kraken2_confidence | number | Kraken2 Confidence score threshold. Default: 0.0. Valid interval: 0-1 | Apply a threshold to determine if a sequence is classified or unclassified. See the [kraken2 manual section on confidence scoring](https://github.com/DerrickWood/kraken2/wiki/Manual#confidence-scoring) for further details about how it works. | 0.0 |


### Minimap2 Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| minimap2filter | string | Filter output of minimap2 by taxids inc. child nodes, E.g. "9606,1404" | Provide a list of taxids if you are only interested in certain ones in your minimap2 analysis outputs. |  |
| minimap2exclude | boolean | Invert minimap2filter and exclude the given taxids instead | Exclude a list of taxids from analysis outputs. | False |
| keep_bam | boolean | Copy bam files into the output directory. |  | False |
| minimap2_by_reference | boolean | Add a table with the mean sequencing depth per reference, standard deviation and coefficient of variation. It adds a scatterplot of the sequencing depth vs. the coverage and a heatmap showing the depth per percentile to the report |  | False |
| min_percent_identity | number | Minimum percentage of identity with the matched reference to define a sequence as classified; sequences with a value lower than this are defined as unclassified. |  | 90 |
| min_ref_coverage | number | Minimum coverage value to define a sequence as classified; sequences with a coverage value lower than this are defined as unclassified. Use this option if you expect reads whose lengths are similar to the references' lengths. |  | 0 |


### Antimicrobial Resistance Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| amr | boolean | Scan reads for antimicrobial resistance or virulence genes | Reads will be scanned using abricate and the chosen database (`--amr_db`) to identify any acquired antimicrobial resistance or virulence genes found present in the dataset. NOTE: It cannot identify mutational resistance genes. | False |
| amr_db | string | Database of antimicrobial resistance or virulence genes to use. |  | resfinder |
| amr_minid | integer | Threshold of required identity to report a match between a gene in the database and fastq reads. Valid interval: 0-100 |  | 80 |
| amr_mincov | integer | Minimum coverage (breadth-of) threshold required to report a match between a gene in the database and fastq reads. Valid interval: 0-100. |  | 80 |


### Report Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| abundance_threshold | number | Remove those taxa whose abundance is equal or lower than the chosen value. | To remove taxa with abundances lower than or equal to a relative value (compared to the total number of reads) use a decimal between 0-1 (1 not inclusive). To remove taxa with abundances lower than or equal to an absolute value, provide a number larger or equal to 1. | 0 |
| n_taxa_barplot | integer | Number of most abundant taxa to be displayed in the barplot. The rest of taxa will be grouped under the "Other" category. |  | 9 |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |
| igv | boolean | Enable IGV visualisation in the EPI2ME Desktop Application by creating the required files. This will cause the workflow to emit the BAM files as well. If using a custom reference, this must be a FASTA file and not a minimap2 MMI format index. |  | False |
| include_read_assignments | boolean | A per sample TSV file that indicates the taxonomy assigned to each sequence. The TSV's will only be output on completion of the workflow and therefore not at all if using the real time option whilst running indefinitely. |  | False |
| output_unclassified | boolean | Output a FASTQ of the unclassified reads. |  | False |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_len | integer | Specify read length lower limit. | Any reads shorter than this limit will not be included in the analysis. | 0 |
| min_read_qual | number | Specify read quality lower limit. | Any reads with a quality lower than this limit will not be included in the analysis. |  |
| max_len | integer | Specify read length upper limit | Any reads longer than this limit will not be included in the analysis. |  |
| threads | integer | Maximum number of CPU threads to use in each parallel workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the number of CPU threads for all such processes. See server threads parameter for Kraken specific threads in the real_time pipeline. | 4 |


