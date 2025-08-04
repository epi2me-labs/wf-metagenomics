#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

include { run_common } from './subworkflows/common_pipeline'
include { run_amr } from './modules/local/amr'
include { minimap_pipeline } from './subworkflows/minimap_pipeline'
// standard kraken2
include { kraken_pipeline } from './subworkflows/kraken_pipeline'

include { prepare_databases } from "./modules/local/databases.nf"
include {
    makeReport;
    getVersions;
    abricateVersion;
} from "./modules/local/common"

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    Pinguscript.ping_start(nextflow, workflow, params)

    // Checking user parameters
    log.info("Checking inputs.")

    // Check maximum and minimum length
    ArrayList fastcat_extra_args = []
    if (params.min_len) { fastcat_extra_args << "-a $params.min_len" }
    if (params.max_len) { fastcat_extra_args << "-b $params.max_len" }
    if (params.min_read_qual) { fastcat_extra_args << "-q $params.min_read_qual" }
    // If BAM files are output, keep runIDs in case they are reused in the wf to track them.
    boolean keep_bam = (params.keep_bam || params.igv)
    if (keep_bam) {fastcat_extra_args << "-H"}

    // Check source param is valid
    sources = params.database_sets
    if (params.containsKey('include_kraken2_assignments')){
        throw new Exception("`include_kraken2_assignments` is now deprecated in favour of `include_read_assignments`.")
    }

    // Stop the pipeline in case not valid parameters combinations
    if (params.classifier == 'minimap2' && params.database) {
        throw new Exception("To use minimap2 with your custom database, you need to use `--reference` (instead of `--database`) and `--ref2taxid`.")
    }

    boolean output_igv = params.igv
    if (params.classifier == 'minimap2' && params.reference && params.igv) {
        ArrayList ref_exts = [".fa", ".fa.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz"]
        if (! ref_exts.any { ext -> file(params.reference).name.endsWith(ext) }) {
            output_igv = false
            log.info("The custom database reference must be a FASTA format file in order to view within IGV.")
        } else {
            output_igv=true
        }
    }

    if (params.classifier == 'kraken2' && params.reference) {
        throw new Exception("To use kraken2 with your custom database, you need to use `--database` (instead of `--reference`) and include the `bracken_dist` within it.")
    }

    // If user provides each database, set to 'custom' the params.database_set
    if (params.reference || params.database) {
        source_name = 'custom'
        // distinguish between taxonomy and database to be able to use taxonomy default db in some cases.
        // this can be potentially risky but might be justified if the reference and ref2taxid use NCBI taxids.
        source_data_database = null
        source_name_taxonomy = params.database_set
        source_data_taxonomy = sources.get(source_name_taxonomy, false)
        log.info("Note: Reference/Database are custom.")
        log.info("Note: Memory available to the workflow must be slightly higher than size of the database $source_name index.")
        if (params.classifier == "kraken2"){
            log.info("Note: Or consider to use the --kraken2_memory_mapping.")
        }

    }
    if(params.taxonomy){
        // this can be useful if the user wants to use a new taxonomy database (maybe updated) but the default reference.
        source_name = params.database_set
        source_data_database = sources.get(source_name, false)
        source_data_taxonomy = null
        log.info("Note: Taxonomy database is custom.")
    } else {
        source_name = params.database_set
        source_data_database = sources.get(source_name, false)
        source_data_taxonomy = sources.get(source_name, false)
        if (!sources.containsKey(source_name) || !source_data_database) {
            keys = sources.keySet()
            throw new Exception("Source $params.database_set is invalid, must be one of $keys")
        }

        def kraken_8GB_dbs = ["Standard-8", "PlusPF-8", "PlusPFP-8"]
        def minimap_dbs = ["ncbi_16s_18s", "ncbi_16s_18s_28s_ITS", "SILVA_138_1"]

        if (!minimap_dbs.contains(source_name) && params.classifier == 'minimap2' && !params.reference){
            throw new Exception("Note: As the classifier parameter is set to `minimap2` the `database_set` parameter must be one of: $minimap_dbs. Use the `kraken2` classifier to access these databases: $kraken_8GB_dbs.")
        }
        if (kraken_8GB_dbs.contains(source_name)){
            log.info("Note: Memory available to the workflow must be slightly higher than size of the database $source_name index (8GB) or consider to use --kraken2_memory_mapping")
        }
    }

    // Input data
    if (params.fastq) {
        ingress_samples = fastq_ingress([
            "input":params.fastq,
            "sample": params.sample,
            "sample_sheet": params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "stats": true,
            "fastcat_extra_args": fastcat_extra_args.join(" "),
            "watch_path": false,
            "per_read_stats": false
        ])
    } else {
        // if we didn't get a `--fastq`, there must have been a `--bam` (as is codified
        // by the schema)
        ingress_samples = xam_ingress([
            "input":params.bam,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "return_fastq": true,
            "keep_unaligned": true,
            "stats": true,
            "watch_path": false,
            "per_read_stats": false
        ])
    }


    // Discard empty samples
    log.info(
        "Note: Empty files or those files whose reads have been discarded after filtering based on " +
        "read length and/or read quality will not appear in the report and will be excluded from subsequent analysis.")
    ingress_samples_filtered = ingress_samples
        | filter { meta, _seqs, _stats ->
            def valid = meta['n_seqs'] > 0
            if (!valid) {
                log.warn "Found empty file for sample '${meta["alias"]}'."
            }
            valid
        }
    // Call the proper pipeline

    // Set minimap2 common options
    ArrayList common_minimap2_opts = [
        "-ax map-ont",
        "--cap-kalloc 100m",
        "--cap-sw-mem 50m",
    ]

    // Run common
    common_versions = getVersions()
    parameters = getParams()

    if (params.exclude_host) {
        host_reference = file(params.exclude_host, checkIfExists: true)
        samples = run_common(ingress_samples_filtered, host_reference, common_minimap2_opts).samples
    } else {
        samples = ingress_samples_filtered
    }

    if (params.classifier == "minimap2") {
        log.info("Minimap2 pipeline.")
        if (keep_bam) {
            common_minimap2_opts = common_minimap2_opts + ["-y"]
        }
        databases = prepare_databases(
            source_data_taxonomy,
            source_data_database
        )
        results = minimap_pipeline(
            samples,
            databases.reference,
            databases.ref2taxid,
            databases.taxonomy,
            databases.taxonomic_rank,
            common_minimap2_opts,
            output_igv
            )
        alignment_stats = results.alignment_reports
    } else {
    // Handle getting kraken2 database files if kraken2 classifier selected
        log.info("Kraken2 pipeline.")
        alignment_stats = Channel.empty()
        databases = prepare_databases(
                source_data_taxonomy,
                source_data_database
        )
        results = kraken_pipeline(
            samples,
            databases.taxonomy,
            databases.database,
            databases.bracken_length,
            databases.taxonomic_rank,
        )
    }

    // Process AMR
    if (params.amr) {
        run_amr = run_amr(
            samples,
            "${params.amr_db}",
            "${params.amr_minid}",
            "${params.amr_mincov}"
        )
        amr_reports = run_amr.reports
        versions = abricateVersion(common_versions)
    } else {
        amr_reports = Channel.empty()
        versions = common_versions
    }
    // Use initial reads stats (after fastcat) QC,
    // and after host_depletion
    // but update meta after running pipelines
    for_report = ingress_samples_filtered
        | map { meta, _path, stats ->
            [ meta.alias, stats ] }
        | combine(
            results.metadata_after_taxonomy,
            by: 0 )  // on alias
        | multiMap { _alias, stats, meta ->
            meta: meta
            stats: stats }
    // Reporting
    makeReport(
        workflow.manifest.version,
        for_report.meta.collect(),
        for_report.stats.collect(),
        results.abundance_table,
        alignment_stats.ifEmpty(OPTIONAL_FILE),
        results.lineages,
        versions,
        parameters,
        databases.taxonomic_rank,
        amr_reports.ifEmpty(OPTIONAL_FILE)
    )
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
