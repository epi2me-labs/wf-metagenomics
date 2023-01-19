# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
- Update kraken databases to latest and ensure relevant taxdump is used.

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

First release.
