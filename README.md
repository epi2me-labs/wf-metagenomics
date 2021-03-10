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

> See the sections below for installation of these prerequisites in various scenarios.
> It is not required to clone or download the git repository in order to run the workflow.

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-metagenomics --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* **report.html** - HTML report document detailing QC metrics and the primary findings of the workflow
* **other.fastq** - fastq of all reads that were classified as something other than human
* **unclassified.fastq** - fastq of all reads that were not classified by centrifuge
* **9606.fastq** - fastq of all reads that classified as human
* **read_classification_master.tsv** - classification result with additional lineage information and fastq stats per read
* **read_classifications.tsv** - classification result for each read - original output of centrifuge


### Supported installations and GridION devices

Installation of the software on a GridION can be performed using the command

`sudo apt install ont-nextflow`

This will install a java runtime, Nextflow and docker. If *docker* has not already been
configured the command below can be used to provide user access to the *docker*
services. Please logout of your computer after this command has been typed.

`sudo usermod -aG docker $USER`

### Installation on Ubuntu devices

For hardware running Ubuntu the following instructions should suffice to install
Nextflow and Docker in order to run the workflow.

1. Install a Jva runtime environment (JRE):

   ```sudo apt install default-jre```

2. Download and install Nextflow may be downloaded from https://www.nextflow.io:

   ```curl -s https://get.nextflow.io | bash```

   This will place a `nextflow` binary in the current working directory, you 
   may wish to move this to a location where it is always accessible, e.g:

   ```sudo mv nextflow /usr/local/bin```

3. Install docker and add the current user to the docker group to enable access:

   ```
   sudo apt install docker.io
   sudo usermod -aG docker $USER
   ```

## Running the workflow

The `wf-metagenomics` workflow can be controlled by the following parameters. The `fastq` parameter
is the most important parameter: it is required to identify the location of the
sequence files to be analysed.

**Parameters:**

- `--reads` specifies path to a FASTQ file (can be tar.gz) (required)
- `--db_path` specifies the directory your centrifuge database directory is in (required)
- `--db_prefix` specifies the name of your centrifuge database (your database files must be in a directory named after your db - see example data) (required)
- `--output-dir` specifies the directory to place your output files in (required)
- `--wfversion` specifies the version of the docker containers to fun when running the workflow in `standard` mode. default: `latest`

To demonstrate the capibilities of the workflow, sample data can be obtained as follows:

1. Sample reads from a mixture of the [Zymo Mock Community]() and a Human cell-line:

   ```
   # Download the sample sequencing reads
   mkdir -p test_data
   wget -O test_data/sample.fastq.gz \
       https://ont-exd-int-s3-euwst1-epi2me-labs.s3-eu-west-1.amazonaws.com/metagenomic_tutorial/sample.fastq.gz
   ```

2. Sample database which will be downloaded and decompressed at `test_data/hpvc/hpvc.tar.gz`.

   ```
   # Download human+viral+prokaryote database
   mkdir -p test_data
   wget -O test_data/hpvc.tar.gz https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/metagenomic_tutorial/hpvc.tar.gz
   # Create a directory and unarchive the files
   cd test_data && mkdir hpvc && cd hpvc && tar -xzvf ../hpvc.tar.gz
   ```

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

```
# run the pipeline with the test data
OUTPUT=output
nextflow run epi2me-labs/wf-metagenomics \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --reads test_data/sample.fastq.gz\
    --db_path test_data/db_store/hpvc \
    --db_prefix hpvc
    --out_dir ${OUTPUT}
    --wfversion latest
```

The output of the pipeline will be found in `./output` for the above
example. This directory contains the nextflow working directories alongside
the primary outputs of the workflow:

#### HTML Report
You can find an example of the report that you will generate here: [link]{} it
gives an overview of the quality of the data along with information about the
classifciations.

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
* label: The name of the output fastq this read was output to e.g. "unclassified" means it was output to `unclassified.fastq`
* len: Length of the read
* meanqscore: Mean qscore of the basecall of the read.


### Running the workflow with Conda

To run the workflow using conda rather than docker, simply replace 

    -profile standard 

with

    -profile conda

in the command above.

### Configuration and tuning

> This section provides some minimal guidance for changing common options, see
> the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for further details.

The default settings for the workflow are described in the configuration file `nextflow.config`
found within the git repository. The default configuration defines an *executor* that will 
use a specified maximum CPU cores (four at the time of writing) and RAM (eight gigabytes).

If the workflow is being run on a device other than a GridION, the available memory and
number of CPUs may be adjusted to the available number of CPU cores. This can be done by
creating a file `my_config.cfg` in the working directory with the following contents:

```
executor {
    $local {
        cpus = 4
        memory = "8 GB"
    }
}
```

and running the workflow providing the `-c` (config) option, e.g.:

```
# run the pipeline with custom configuration
nextflow run epi2me-labs/wf-metagenomics \
    -c my_config.cfg \
    ...
```

The contents of the `my_config.cfg` file will override the contents of the default
configuration file. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)
for more information concerning customized configuration.

**Using a fixed conda environment**

By default, Nextflow will attempt to create a fresh conda environment for any new
analysis (for reasons of reproducibility). This may be undesirable if many analyses
are being run. To avoid the situation a fixed conda environment can be used for all
analyses by creating a custom config with the following stanza:

```
profiles {
    // profile using conda environments rather than docker
    // containers
    fixed_conda {
        docker {
            enabled = false
        }
        process {
            withLabel:artic {
                conda = "/path/to/my/conda/environment"
            }
            shell = ['/bin/bash', '-euo', 'pipefail']
        }
    }
}
```

and running nextflow by setting the profile to `fixed_conda`:

```
nextflow run epi2me-labs/wf-metagenomics \
    -c my_config.cfg \
    -profile fixed_conda \
    ...
```

## Updating the workflow

Periodically when running the workflow, users may find that a message is displayed
indicating that an update to the workflow is available.

To update the workflow simply run:

    nextflow pull epi2me-labs/wf-metagenomics

## Building the docker container from source

The docker image used for running the `wf-metagenomics` workflow is available on
[dockerhub](https://hub.docker.com/repository/docker/ontresearch/wf-metagenomics).
The image is built from the Dockerfile present in the git repository. Users
wishing to modify and build the image can do so with:

```
CONTAINER_TAG=ontresearch/wf-metagenomics:latest

git clone https://github.com/epi2me-labs/wf-metagenomics
cd wf-metagenomics

docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=ontresearch/base-workflow-image:v0.1.0 \
    .
```

In order to run the workflow with this new image it is required to give
`nextflow` the `--wfversion` parameter:

```
nextflow run epi2me-labs/wf-metagenomics \
    --wfversion latest
```

## Useful links

* [centrifuge](https://ccb.jhu.edu/software/centrifuge/)
* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
