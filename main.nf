#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
include { minimap_pipeline } from './subworkflows/minimap_pipeline'
include { kraken_pipeline } from './subworkflows/kraken_pipeline'

nextflow.preview.recursion=true
   
// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfmetagenomics"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
       
    }

    if (params.kraken2 && params.minimap2){
        log.info("")
        log.info("You must specify a classification method(s) with --kraken2 or --minimap2 but not both")
        exit 1
    }

    dataDir = projectDir + '/data'

    // Ready the optional file
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    // Acquire report template
    template = file("$projectDir/bin/report-visualisation.html")

    // Checking user parameters
    log.info("Checking inputs.")

    if (!(params.minimap2 || params.kraken2)) {
        log.info("")
        log.info("You must specify a classification method with --kraken2 or --minimap2")
        exit 1
    }

    if (params.kraken2filter && !params.kraken2) {
        log.info("")
        log.info("Usage of kraken2filter requires `--kraken2`.")
        exit 1
    }

    // Check source param is valid
    sources = params.sources
    source_name = params.source
    source_data = sources.get(source_name, false)
    if (!sources.containsKey(source_name) || !source_data) {
        keys = sources.keySet()
        log.info("Source $params.source is invalid, must be one of $keys")
        exit 1
    }

    // Grab taxonomy files
    taxonomy = file(sources[source_name]["taxonomy"], type: "file")
    if (params.taxonomy) {
        log.info("Checking custom taxonomy mapping exists")
        taxonomy = file(params.taxonomy, type: "dir", checkIfExists:true)
    }

    // Handle getting alignment reference files if minimap2 is enabled
    reference = null
    refindex  = null
    ref2taxid = null
    if (params.minimap2) {
        // .fasta
        if (params.reference) {
            log.info("Checking custom reference exists")
            reference = file(params.reference, type: "file", checkIfExists:true)
        } else {
            source_reference = source_data.get("reference", false)
            if (!source_reference) {
                log.info(
                    "Error: Source $source_name does not include a reference for "
                    + "use with minimap2, please choose another source, "
                    + "provide a custom reference or disable minimap2.")
                exit 1
            }
            reference = file(source_reference, type: "file")
        }
        // .fasta.fai
        refindex = file(sources[source_name]["refindex"], type: "file")
        if (params.reference) {
            log.info("Checking custom reference index exists")
            refindex = file(params.reference + '.fai', type: "file")
            if (!refindex.exists()) {
                refindex = file(OPTIONAL, type: "file")
            }
        }
        // .ref2taxid.csv
        ref2taxid = file(sources[source_name]["ref2taxid"], type: "file")
        if (params.ref2taxid) {
            log.info("Checking custom ref2taxid mapping exists")
            ref2taxid = file(params.ref2taxid, type: "file", checkIfExists:true)
        }

        samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet])


        results = minimap_pipeline(
        samples, reference, refindex, ref2taxid, taxonomy,
        template)

    }

    // Handle getting kraken2 database files if kraken2 is enabled
    database = null
    kmer_distribution = null
    if (params.kraken2) {
        // kraken2.tar.gz
        if (params.database) {
            log.info("Checking custom kraken2 database exists")
            database = file(params.database, type: "dir", checkIfExists:true)
        } else {
            source_database = source_data.get("database", false)
            if (!source_database) {
                log.info(
                    "Error: Source $source_name does not include a database for "
                    + "use with kraken2, please choose another source, "
                    + "provide a custom database or disable kraken2.")
                exit 1
            }
            database = file(source_database, type: "file")
        }
        // kmer_distrib
        kmer_dist_path = sources[source_name]["kmer_dist"]
        if (params.bracken_dist) {
            log.info("Checking custom bracken2 database exists")
            kmer_distribution = file(
                params.bracken_dist, type: "file", checkIfExists:true)
        } else if (!kmer_dist_path) {
            kmer_distribution = file(OPTIONAL, type: "file")
            if (params.kraken2bracken) {
                log.info("Kmer distribution not found, bracken2 disabled.")
            }
        } else {
            kmer_distribution = file(sources[source_name]["kmer_dist"], type: "file")
        }
        // check combination of params are set
        if (params.watch_path && !params.run_indefinitely){
            if (!params.read_limit){
                log.error("Must specify read limit parameter if run indefinitely is false")
                exit 1 
            }
            log.info("Workflow will stop processing files after ${params.read_limit} reads when run_indefinitely is set to False")  
        }
        results = kraken_pipeline(
        reference, refindex, ref2taxid, taxonomy,
        database, kmer_distribution, template)

    }

    output(results)
}


if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }
    
    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
