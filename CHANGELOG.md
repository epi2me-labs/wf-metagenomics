# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [Unreleased]
### Fixed
- Switch to markdown links in the outputs table in the README.
- Exclude samples if all the reads are removed during host depletion.
### Added
- `igv` option to enable IGV in the EPI2ME Desktop Application.
- `include_read_assignments` option to output a file with the taxonomy of each read.
- `Reads` section in the report to track the number of reads after filtering, host depletion and unclassified.
### Changed
- `keep_bam` is now only required to output BAM files.
- `include_kraken2_assignments` has been replaced by `include_read_assignments`.
- Update databases:
    - Kraken2 indexes to the ones released 2024-06-05
    - Taxonomy database to the one released 2024-09-01
### Removed
- `split-prefix` parameter, as the workflow automatically enables this option for large reference genomes.

## [v2.10.1]
### Fixed
- Fix the json.decoder.JSONDecodeError when running with `--amr`.

## [v2.10.0]
### Changed
- The workflow now uses the `fastcat` read length and quality histograms instead of the per-read stats in the report process.
### Added
- Output IGV configuration file if the `keep_bam` option is enabled and a custom reference is provided (in minimap2 mode).
- Output reduced reference file if the `keep_bam` option is enabled (in minimap2 mode).
- `abundance_threshold` reduces the number of references to be displayed in IGV.
### Fixed
- Complete previous taxonomic rank when there are no parent nodes.
- Checking the correspondence between the reference and ref2taxid now also works with compressed references.
- `exclude-host` can input a file in the EPI2ME Desktop Application.

## [v2.9.4]
### Added
- Statistics derived from minimap2 alignment are now in the output if the `keep_bam` option is enabled.
- Reads below percentages of identity (`min_percent_identity`) and the reference covered (`min_ref_coverage`) are considered as unclassified in the minimap2 approach.
### Fixed
- Request less memory if `kraken2_memory_mapping` is used.
- Show the percentage of each species when hovering over the taxonomy bar plot.
### Changed
- `bam` folder within output has been renamed to `bams`.

## [v2.9.3]
### Fixed
- Files that are empty following the fastcat filtering are discarded from downstream analyses.

## [v2.9.2]
### Fixed
- "Can only use .dt accessor with datetimelike values" error in makeReport 
- "invalid literal for int() with base 10" error in makeReport

## [v2.9.1]
### Fixed
- Avoid argument list lengths that may be too long when using glob patterns.
- The Kraken2 pipeline sometimes reporting wrong numbers for unclassified counts.

## [v2.9.0]
### Fixed
- Minimum memory requirement for custom databases depends on the database size.
### Added
- Workflow now accepts BAM or FASTQ files as input (using the `--bam` or `--fastq` parameters, respectively).
- Run ID's now included in the output BAM files.
- Add nextflow tags to keep track of the samples along the different processes.
- Output BAM files without host reads.
### Changed
- Default for `--n_taxa_barplot` increased from 8 to 9.
- Replace the eCharts barplot with an equivalent Bokeh version.
- `--database_set` parameter is now `Standard-8` when the `--classifier` parameter is set to `kraken2`.

## [v2.8.0]
### Added
- Standard 8GB kraken2 database.
- Update docs.
### Fixed
- Heatmap generated when `--minimap2_by_reference` is enabled references with a mean coverage of less than 1% of the one with the largest value are omitted.

## [v2.7.0]
### Fixed
- Use store_dir without staging files from the web. Kraken2 can run offline if the databases have been previously stored.
- Fastcat plots showing the stats in the report before removing host sequences when `--exclude_host` in the minimap2 pipeline.
- Real time kraken workflow hanging indefinitely when attempting to start kraken server with too many threads.

### Added
- Minimap pipeline is also able to use store_dir to store databases and run offline if the databases have been previously stored.
- Kraken2 pipeline accepts a sample sheet if the real time option is disabled.
- Only taxa present in the abundance table above the `--abundance_threshold` will appear in the alignment summary table (which is only generated when `--minimap2_by_reference` is enabled). 

### Removed
- `--bracken_dist`: the bracken additional file for the database must be included in the database folder, as it is in the kraken2 indexes and when the database is generated.
- Default local executor CPU and RAM limits

### Changed
- `--watch_path` is now called `--real_time` and enables the kraken2 pipeline to classify reads as they are written with watch_path.
- The kraken2 workflow can now be used without `--real_time`, this will use the serverless kraken2 executable.
- Barcode directories must now be named in the format `barcodeNN`, where NN is at least two digits (e.g. `barcode01`).
- Barcode directories must now have the same number of characters (e.g. `barcode01` cannot be provided with `barcode001`).

## [v2.6.1]
### Fixed
- Broken report when the dataframe is filtered using the `--abundance_threshold`.
- Taxonomy abundances barplot was not showing more abundant species.

## [v2.6.0]
### Fixed
- Broken plots caused by single quotes in NCBI taxon names.
### Added
- Add the abundance_table_rank.tsv in the output for the last analysed rank.
- Optional `--minimap2_by_reference` parameter to output the sequencing depth and coverage of each matched reference in the database.

## [v2.5.0]
### Added
- `--kraken_confidence` to specify a threshold score.
- `--exclude_host`: Optional parameter can accept a FASTA/MMI file with a host reference to be excluded from the analysis.
- `--include_kraken2_assignments`: Output the classification of each read.

## [v2.4.4]
### Added
- `--abundance_threshold`: filter taxa based on their abundances.
- `--n_taxa_barplot`: control the number of taxa displayed in the barplot.
- Plot the taxa abundance distribution (e.g. Species abundance distribution plots).

## [v2.4.3]
### Changed
- Remove abricate version if AMR does not run.
### Fixed
- Changelog format.
### Added
- Alpha diversity indices: Berger-Parker dominance index, Fisher’s alpha.

## [v2.4.2]
### Changed
-Bumped minimum required Nextflow version to 23.04.2.

## [v2.4.1]
### Fixed
- Kraken2 pipeline: all the samples are shown in the report.
### Changed
- Any sample aliases that contain spaces will be replaced with underscores.

## [v2.4.0]
### Added
- Antimicrobial resistance gene identification using Abricate.

## [v2.3.0]
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
- A new option `--min_read_qual` to filter by quality score.
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

