# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Changed
- Better help text on cli

## [v1.1.4]

## Added
- pluspf8, ncbi_16s_18s_28s_ITS databases
- Add all sample tool combinations to report
## Changed
- Enable kraken2 by default
- Clarify error messages
## Fixed
- Handle no assignments bracken error

## [v1.1.3]

## Added
- New docs format.
- Render bokeh.

## [v1.1.2]

## Changed
- Update nextflow_schema.json

## [v1.1.1]

## Fixed
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
