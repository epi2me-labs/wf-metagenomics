# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- A new option `kraken2_memory_mapping` to avoid kraken2 loading the database into process-local RAM.
- `--keep_bam` parameter to write BAM files into the output directory (minimap pipeline).
- Lineages sunburst plot added to the report.
- SILVA.138 database available for both kraken2 and minimap2 pipelines.

### Changed
- `bracken_level` parameter has been replaced by `taxonomic_rank` to choose the taxonomic rank at which to perform the analysis. It works in both pipelines.
- Updated example command displayed when running `--help`.
- Updated GitHub issue templates to force capture of more information.
- Bumped minimum required Nextflow version to 22.10.8.
- Enum choices are enumerated in the `--help` output.
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice.

### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`.
- Extract reads using `--minimap2filter` and `--minimap2exclude` filters. The extracted reads are in the output/filtered folder.

## [v2.2.1]
### Added
- A new option `min_read_qual` to filter by quality score.
- Configuration for running demo data in AWS
- AWS configuration for external kraken2_server for demonstration at LC23

### Fixed
- Fix minimum and maximum length read filter.

### Removed
- Default region and AWS CLI path for AWS batch profile

## [v2.2.0]
### Changed
- Updated existing databases.
- Docker will use an ARM platform image on appropriate devices.

### Added
- A new PFP-8 database option.
- New test_data with Bacteria, Archaea and Fungi.

### Fixed
- Fix file names when exporting tables.

## [v2.1.1]
### Fixed
- Include 'kingdom' for Eukarya.
- Add ability to use an external kraken2 server.

## [v2.1.0]
### Updated
- New fastqingress.
- New report with ezcharts.

### Added
- Stacked barplot for most abundant taxa.
- Show rank information in abundance tables.
- Export function from tables.

### Fixed
- Fix crash in the report with one sequence in the fastq.
- Use kraken2 with parallelization in single client.

## [v2.0.10]
### Fixed
- Remove symbolic links from store_dir.

## [v2.0.9]
### Changed
- Update kraken databases to latest and ensure relevant taxdump is used.

### Added
- Plot species richness curves.
- Provide (original and rarefied, i.e. all the samples have the same number of reads) abundance tables listing taxa per sample for a given taxonomic rank.
- Add diversity indices.
- Memory requirement help text.

## [v2.0.8]
### Fixed
- Example_cmd in config.
- Minimap2 subworkflow fixed for when no alignments.
- Remove quality 10 parameter from Minimap2 subworkflow.

## [v2.0.7]
### Fixed
- Issue where processing more than ~28 input files lead to excessive memory use.
- Minor typos in docs.

## [v2.0.6]
### Changed
- Updated description in manifest

## [v2.0.5]
### Fixed
- The version in the config manifest is up to date.

## [v2.0.4]
### Fixed
- Issue where discrepancies between taxonomy and databases led to error.

### Added
- `nextflow run epi2me-labs/wf-metagenomics --version` will now print the workflow version number and exit.

### Changed
- Parameter name for selecting known database is now `--database_set` (was `--source`).
- Add classifier parameter and only allow running of minimap2 or kraken2 workflow.
- Workflow logic in kraken workflow has been reorganised for simpler parallelism.

### Removed
- `-profile conda` is no longer supported, users should use `-profile standard` (Docker) or `-profile singularity` instead
- `--run_indefinitely` parameter removed, instead implied when `--read_limit` set to null.

## [v2.0.3]
### Fixed
- Add a test and fix for if all files in one directory are unclassified
- Check if fastq input exists

## [v2.0.2]
### Changed
- Use store directory for database.
- Use per file kraken_report instead of cumulative.
- Kraken2-server v0.0.8.

### Added
- Add a run indefinitely parameter.

### Fixed
- Batch size breaking fastcat step.
- Consider white space in bracken report.
- Handling for unclassified with Bracken.

## [v2.0.1]
### Fixed
- Handling with kraken2 for single input file

### Removed
- Removed sanitize option

## [v2.0.0]
### Fixed
- Output argument in Fastqingress homogenised.

### Changed
- Bumped base container to v0.2.0
- Kraken workflow now in real time mode with watch_path
- Kraken and Minimap now in subworkflows
- Fastqingress metadata map
- Can only run Kraken or Minimap subworkflow not both
- Better help text on cli
- Fastq ingress and Args update
- Set out_dir option type to ensure output is written to correct directory on Windows

## [v1.1.4]
### Added
- pluspf8, ncbi_16s_18s_28s_ITS databases
- Add all sample tool combinations to report

### Changed
- Enable kraken2 by default
- Clarify error messages

### Fixed
- Handle no assignments bracken error

## [v1.1.3]
### Added
- New docs format.
- Render bokeh.

## [v1.1.2]
### Changed
- Update nextflow_schema.json

## [v1.1.1]
### Fixed
- Overriding taxonomy now works correctly
- Added missing threads param to kraken2

## [v1.1.0]
### Added
- Report now includes dynamic sankey visualisation and table
- Nextflow schema

### Changed
- Updated to use new fastqingress module, permitting single .fastq input
- Rewired DAG to rely on sample id's rather than filenames

### Fixed
- Handle bracken failure when there are no classifications
- Handle cyclic dag issue when taxonomy has duplicate names

## [v1.0.0]
* First release.

