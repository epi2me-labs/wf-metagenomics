# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- Help message when required parameters are not given.
- `--split` to customise partitioning of output fastq data

### Changed
- Now uses `ncbitaxonomy` rather than `ete3`
- Outputs dirs within `fastq_bundles` containing reads partitioned by 
  rules specified by `--split`

## [v0.1.0]
### Added
- Help message when required parameters are not given
- `--download` command to download other sample metagenomic databases
- `--threads` to change number of threads given to centrifuge

### Changed
- README example command now refers to the zymo db rather hvpc
- `--reads` param changed to `--fastq`
- Performance improvement to `split_by_master.py`

## [v0.0.3]

### Removed
- Superfluous "about" report section 

## [v0.0.2]

### Added
- New minimal centrifuge database included with the project
- New minimal sample dataset for testing

### Changed
- Report reads stats section using standard aplanat report section

### Removed
- Report section comparing unclassified vs classified reads removed

## [v0.0.1]

Initial release

### Added
- Basic running of centrifuge, manipulation of outputs, and reporting.
