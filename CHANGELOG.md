# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
