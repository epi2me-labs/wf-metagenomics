# Metagenomic classifications

This repository contains a [nextflow](https://www.nextflow.io/) workflow
performing metagenomic classification of Oxford Nanopore Technologies'
sequencing datasets. The workflow will perform the classification and
produce useful groupings of the reads for further analysis.

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

* **report.html** - HTML report document detailing QC metrics and the primary findings of the workflow
* **fastq_budles** - Directory containing fastqs grouped by `--split` opt. By default this 
  will be a directories: `fungi`, `bacteria`, `viral` containing fastqs grouped by phylum
  and further `unclassified` containining reads that were not classified and `else` which
  contains, split by superkingdom, all other read classifications.
* **seqs.txt** - output of `fastcat` containing per read stats
* **seqs.fastq** - concatenated fastq of all reads passed to workflow 
* **read_classifications.tsv** - classification result for each read - original output of centrifuge
* **read_classification_master.tsv** - classification result with additional lineage information and fastq stats per read

## Running the workflow

The `wf-metagenomics` workflow can be controlled by the following parameters.

**Parameters:**

- `--fastq` specifies path to a FASTQ file (can be tar.gz) (required) e.g. `test_data/sample.fastq.gz` or `test_data`
- `--db_path` specifies the directory your centrifuge database files (*.cf) are in (required) e.g. `test_data/db_store/zymo/`
- `--db_prefix` specifies the name of your centrifuge database (required) `zymo` for `zymo.*.cf` database
- `--out_dir` specifies the directory to place your output files in (required). default: `output`
- `--wfversion` specifies the version of the docker containers to fun when running the workflow in `standard` mode. default: `latest`
- `--threads` specifies the number of threads available to centrifuge (default: 8)
- `--split` add "rules" (see section below for details) for labelling / splitting results (default: splitting by superkingdom)

To demonstrate the capabilities of the workflow sample data has been included.  A selection of sample reads from a mixture 
of the [Zymo Mock Community](https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards/products/zymobiomics-microbial-community-dna-standard) 
and a Human cell-line are included along with a sample database containing **only** 
Zymo mock community references is included at `/test_data/db_store/zymo.tar.gz`.

Before running the workflow please decompress the Zymo centrifuge database:
   ```
   # Decompress the Zymo centrifuge database
   tar xvzf test_data/db_store/zymo.tar.gz -C test_data/db_store/
   ```
   
> Please note that the example database is only to act as a very minimal example database
> and is not suitable for analysis.  Please follow the instructions at the bottom of this README
> to download a more comprehensive database containing a wider range of organisms.

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run` to run with conda use `-profile conda`:

```
# run the pipeline with the test data
OUTPUT=output
nextflow run epi2me-labs/wf-metagenomics \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq test_data/sample.fastq.gz\
    --db_path test_data/db_store/zymo \
    --db_prefix zymo
    --out_dir ${OUTPUT}
    --wfversion latest
```

The output of the pipeline will be found in `./output` for the above
example. This directory contains the nextflow working directories alongside
the primary outputs of the workflow:

#### .fastq files

Three fastq files are created: `9606.fastq`, `other.fastq` and
`unclassified.fastq`. These three files contain all the reads that were
processed by the workflow but split by the classification that centrifuge gave.
`9606.fastq` contains reads that were classified as human (9606 is the taxid of
the genus Homo). Reads that classified as a particular organism but not Human
(or at a higher rank than genus within the same lineage) are output into
`other.fastq` and `unclassified.fastq` contains all reads that were not
classified.  Although most experiments will output reads that cannot be
classified either because of preparaion/quality issues or because the organisms
that were present were not sufficiently represented in the database given.
Please note, no classifications can be made for references that are not present
in the database as it compares against what it "knows".

#### `read_classification_master`

This file details for each read the following fields. For up-to-date
documentation on the centrifuge output fields please see the centrifuge
documentation:
[here](https://ccb.jhu.edu/software/centrifuge/manual.shtml#centrifuge-classification-output)

* readID - The unique identifyer for the read present in the fastq file
* centrifuge output:
    > seqID, 
    > taxID, 
    > score, 
    > 2ndBestScore
    > hitLength
    > queryLength
    > numMatches
* lineage: superkingdom, kingdom, phylum, class, order, family, genus, species,
* label: The name of the file this fastq was output to e.g. "unclassified" means it was output to `unclassified.fastq`
* len: Length of the read
* meanqscore: Mean qscore of the basecall of the read.

### Rules

The `--split` parameter can be used to label and split results into  particular clades.
The syntax is 1 or more space separated rules.

```
    --split "label:action:taxid label:action:taxid"
```

Each rule is applied on the remaining reads that are not already sorted so it's 
important to use more specific rules before more generalised ones.

Individual rules are colon-separated as label:action:taxid. 

* `label` indicates the folder the file that the individual fastq record will be stored in.
* `action` indicates the way the fastqs within a particular labelled folder should be grouped. 
* `taxid` refers to the taxid and those below that should be included within 
this action.

#### Actions
`action` can take one of a specific list of actions detailed below:

##### Group by clade:

* `superkingdom` - split by superkingdom of the classification
* `kingdom` - split by kingdom of the classification
* `phylum` - split by phylum of the classification 
* `class` - split by class of the classification 
* `order` - split by order of the classification 
* `family` - split by family of the classification 
* `genus` - split by genus of the classification 
* `species` - split by species of the classification 

If a read is classified at a level below the grouping `action` selected e.g. it's 
classified at the species level but the grouping `action` is at the class level:  
the fastq record will be stored in a file named after the `class` within the 
classification that read's classified lineage.  However if a classified taxid is 
above the classification rank e.g. the read is classified as "Bacteria" but the 
action `class` is used then this read will be placed in `{label}/other.fastq`.

##### Group by taxid:

* `taxid` - split by taxid classification

The read will be grouped together with other reads that are classified as exactly
the same taxid e.g. *E.coli* reads classified at the species level (taxid=562) will 
be grouped together in one file (562.fastq) and *E.coli* reads classified at the 
strain level e.g. *E. coli O17 str. K12a* (taxid=1010810) will be grouped together in another file
e.g. 1010810.fastq within the same folder.

#### Examples

##### All bacteria in a "Contamination" folder

If you wanted all the reads classified as Bacteria (taxid=2)
to be included in a folder called "Contamination" and that these 
fastq should be grouped at the "family" such that you would have a separate file 
for fastqs classified as "Enterobacteriaceae" (taxid=543) and 
"Lactobacillaceae" (taxid=33958) you would get the following:
```
# redacted nextflow command with split
nextflow [...] --split "Contamination:family:2"
```
Results in the following files:
```
Contamination/543.fastq # Enterobacteriaceae and below 
Contamination/33958.fastq # Lactobacillaceae and below
Contamination/other.fastq # Reads classified as bacteria but above the family level
other.fastq # Reads classified without Bacteria (taxid=2) within their lineage
unclassified.fastq # Reads that were not classified
```

##### All classifications in one "Classified" folder
```
# redacted nextflow command with split
nextflow [...] --split "Classified:superkingdon:1"
nextflow [...] --split "Classified:superkingdon:"
```
If you want all reads in a folder called "Classified" you can use 1 as the 
root taxid or leave it blank. This will result in the following folders:

```
Classified/2.fastq # Bacteria
Classified/2759.fastq # Eukaryota
unclassified.fastq # unclassified reads
```

##### All E.coli in one folder and everything else in Classified folder
```
# redacted nextflow command with split
nextflow [...] --split "Escherichia:taxid:561 Classified:superkingdon:1"
nextflow [...] --split "Escherichia:taxid:561 Classified:superkingdon:"
```
Files created:
```
Ecoli/562.fastq # Escherichia coli
Ecoli/1010810.fastq # Escherichia coli O17 str. K12a
Ecoli/1355083.fastq # Escherichia coli ZH193
Classified/2.fastq # Bacteria
Classified/2759.fastq # Eukaryota
unclassified.fastq
```

## Note on databases

Note: This database only contains the Zymo mock community references so is only useful as a lightweight example.
You can download a more comprehensive database containing refseq human, viral and prokaryotic references (`hvpc`) using 
the instructions below. Please bear in mind that this reference database is larger (27MB vs 20GB compressed) so may take some time to download 
and will require a significant amount of memory (36GB RAM) to run successfully.

The `hpvc` sample database which will be downloaded and decompressed at `test_data/hpvc.*.cf`.

   ```
   # Download human+viral+prokaryote+covid database
   nextflow run main.nf -w output/workspace -profile standard --db_path test_data/db_store --db_prefix hpvc --wfversion latest --download
   # Run the sample dataset through the expanded database
   nextflow run main.nf -w output/workspace -profile standard --fastq test_data/sample.fastq.gz --db_path test_data/db_store --db_prefix hpvc --out_dir output --wfversion latest
   ```

## Useful links

* [centrifuge](https://ccb.jhu.edu/software/centrifuge/)
* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
