#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 
 
 
process unpackDatabase {
    label "wfmetagenomics"
    cpus 1
    input:
        path database
        path kmer_distribution
    output:
        path "database_dir"
    """
    if [[ $database == *.tar.gz ]]
    then
        mkdir database_dir
        tar xf "${database}" -C database_dir
    elif [ -d "${database}" ]
    then
        mv "${database}" database_dir
    else
        echo "Error: database is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    mv "${kmer_distribution}" database_dir
    """
}


process unpackTaxonomy {
    label "wfmetagenomics"
    cpus 1
    input:
        path taxonomy
    output:
        path "taxonomy_dir"
    """
    if [[ "${taxonomy}" == *.tar.gz ]]
    then
        mkdir taxonomy_dir
        tar xf "${taxonomy}" -C taxonomy_dir
    elif [ -d "${taxonomy}" ]
    then
        mv "${taxonomy}" taxonomy_dir
    else
        echo "Error: taxonomy is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}


process combineFilterFastq {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        tuple(
            val(sample_id),
            path("${sample_id}.fastq"),
            emit: filtered)
        tuple(
            val(sample_id),
            path("${sample_id}.stats"),
            emit: stats)
    shell:
    """
    fastcat \
        -a "${params.min_len}" \
        -b "${params.max_len}" \
        -q 10 \
        -s "${sample_id}" \
        -r "${sample_id}.stats" \
        -x "${directory}" > "${sample_id}.fastq"
    """
}


process minimap2 {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path(reads)
        path reference
        path refindex
        path ref2taxid
        path taxonomy
    output:
        tuple(
            val(sample_id),
            path("*.bam"),
            path("*.bam.bai"),
            emit: bam)
        tuple(
            val(sample_id),
            path("*assignments.tsv"),
            emit: assignments)
        tuple(
            val(sample_id),
            path("*lineages.txt"),
            emit: lineage_txt)
        tuple(
            val(sample_id),
            path("*lineages.json"),
            emit: lineage_json)
    script:
        def split = params.split_prefix ? '--split-prefix tmp' : ''
    """
    minimap2 -t "${params.threads}" -ax map-ont "${split}" "${reference}" "${reads}" \
    | samtools view -h -F 2304 - \
    | format_minimap2.py - -o "${sample_id}.minimap2.assignments.tsv" -r "${ref2taxid}" \
    | samtools sort -o "${sample_id}.bam" -
    samtools index "${sample_id}.bam"
    awk -F '\\t' '{print \$3}' "${sample_id}.minimap2.assignments.tsv" > taxids.tmp
    taxonkit \
        --data-dir "${taxonomy}" \
        lineage -R taxids.tmp \
        | aggregate_lineages.py -p "${sample_id}.minimap2"
    """
}


process extractMinimap2Reads {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path(bam), file(bai)
        path ref2taxid
        path taxonomy
    output:
        tuple(
            val(sample_id),
            path("*extracted.fastq"), 
            emit: extracted)
    script:
        def policy = params.minimap2exclude ? '--exclude' : ''
    """
    taxonkit \
        --data-dir "${taxonomy}" \
        list -i "${params.minimap2filter}" \
        --indent "" > taxids.tmp
    extract_minimap2_reads.py \
        "${bam}" \
        -r "${ref2taxid}" \
        -o "${sample_id}.minimap2.extracted.fastq" \
        -t taxids.tmp \
        "${policy}"
    """
}


process kraken2 {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path(reads)
        path database
        path taxonomy
    output:
        tuple(
            val(sample_id),
            path("*.classified.fastq"),
            emit: classified)
        tuple(
            val(sample_id),
            path("*.unclassified.fastq"),
            emit: unclassified)
        tuple(
            val(sample_id),
            path("*assignments.tsv"),
            emit: assignments)
        tuple(
            val(sample_id),
            path("*lineages.txt"),
            emit: lineage_txt)
        tuple(
            val(sample_id),
            path("*lineages.json"),
            emit: lineage_json)
        tuple(
            val(sample_id),
            path("*kraken2_report.txt"),
            emit: kraken2_report)
    """
    kraken2 \
        --db "${database}" \
        --threads "${params.threads}" \
        --report "${sample_id}.kraken2_report.txt" \
        --classified-out "${sample_id}.kraken2.classified.fastq" \
        --unclassified-out "${sample_id}.kraken2.unclassified.fastq" \
        "${reads}" > "${sample_id}.kraken2.assignments.tsv"
    awk -F '\\t' '{print \$3}' "${sample_id}.kraken2.assignments.tsv" > taxids.tmp
    taxonkit \
        --data-dir "${taxonomy}" \
        lineage -R taxids.tmp \
        | aggregate_lineages.py -p "${sample_id}.kraken2"
    """
}


process extractKraken2Reads {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path(reads), path(kraken_assignments), path(kraken_report)
    output:
        tuple(
            val(sample_id),
            path("*extracted.fastq"), 
            emit: extracted)
    script:
        def taxids = (params.kraken2filter as String).replaceAll(',',' ')
        def policy = params.kraken2exclude ? '--exclude' : ''
    """
    extract_kraken_reads.py \
        -k "${kraken_assignments}" \
        -r "${kraken_report}" \
        -s1 "${reads}" \
        -o "${sample_id}.kraken2.extracted.fastq" \
        -t "${taxids}" \
        --fastq-output \
        --include-children
        "${policy}"
    """
}


process bracken {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path(kraken2_report)
        path database
    output:
        tuple(
            val(sample_id),
            path("*bracken_report.txt"), 
            emit: bracken_report)
    """
    run_bracken.py \
        "${database}" \
        "${kraken2_report}" \
        "${params.bracken_length}" \
        "${params.bracken_level}" \
        "${sample_id}.bracken_report.txt"
    """
}


process getVersions {
    label "wfmetagenomics"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    taxonkit version | sed 's/ /,/' >> versions.txt
    kraken2 --version | head -n 1 | sed 's/ version /,/' >> versions.txt
    """
}


process getParams {
    label "wfmetagenomics"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process makeReport {
    label "wfmetagenomics"
    input:
        path stats
        path lineages
        path "versions/*"
        path "params.json"
        path template
    output:
        path "wf-metagenomics-*.html"
    script:
        report_name = "wf-metagenomics-" + params.report_name + '.html'
    """
    report.py \
        "${report_name}" \
        --versions versions \
        --params params.json \
        --summaries "${stats}" \
        --lineages "${lineages}" \
        --vistempl "${template}"
    """
}


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


// workflow module
workflow pipeline {
    take:
        samples
        reference
        refindex
        ref2taxid
        taxonomy
        database
        kmer_distribution
        template
    main:
        outputs = []
        lineages = Channel.empty()
        taxonomy = unpackTaxonomy(taxonomy)
        if (params.kraken2) {
            database = unpackDatabase(
                database,
                kmer_distribution
            )
        }

        // Initial reads QC
        reads = combineFilterFastq(samples)

        // Run Kraken2
        if (params.kraken2) {
            kr2 = kraken2(
                reads.filtered, 
                database, 
                taxonomy
            )
            reads_to_align = kr2.classified
            lineages = lineages.mix(kr2.lineage_json)
            outputs += [
                kr2.classified,
                kr2.unclassified,
                kr2.lineage_txt,
                kr2.assignments,
                kr2.kraken2_report
            ]
            if (params.kraken2bracken
                && kmer_distribution.name != 'OPTIONAL_FILE') {
                br = bracken(
                    kr2.kraken2_report,
                    database
                )
                outputs += [br.bracken_report]
            }
            if (params.kraken2filter) {
                kr2_filt = extractKraken2Reads(
                    reads_to_align,
                    kr2.assignments,
                    kr2.kraken2_report
                )
                outputs += [kr2_filt.extracted]
                reads_to_align = kr2_filt.extracted
            }
        } else {
            reads_to_align = reads.filtered
        }

        // Run Minimap2
        if (params.minimap2) {
            mm2 = minimap2(
                reads_to_align,
                reference,
                refindex,
                ref2taxid,
                taxonomy
            )
            lineages = lineages.mix(mm2.lineage_json)
            outputs += [
                mm2.bam,
                mm2.assignments,
                mm2.lineage_txt
            ]
            if (params.minimap2filter) {
                mm2_filt = extractMinimap2Reads(
                    mm2.bam,
                    ref2taxid,
                    taxonomy
                )
                outputs += [mm2_filt.extracted]
            }
        }

        // Reporting
        software_versions = getVersions()
        workflow_params = getParams()
        report = makeReport(
            reads.stats.flatMap { it -> [ it[1] ] }.collect(),
            lineages.flatMap { it -> [ it[1] ] }.collect(),
            software_versions.collect(), 
            workflow_params,
            template
        )
    emit:
        report.concat(
            *outputs.collect { chan ->
                chan.flatMap { it -> [ it[1] ] }},
            software_versions,
            workflow_params
        )
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)

workflow {
    dataDir = projectDir + '/data'

    // Ready the optional file
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    // Acquire report template
    template = file("$projectDir/bin/report-visualisation.html")

    // Checking user parameters
    log.info("Checking inputs.")

    if (!(params.minimap2 || params.kraken2)) {
        log.info("")
        log.info("You must specify a classification method(s) with --kraken2 or --minimap2 (or both)")
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
    }

    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet,
        params.sanitize_fastq)

    results = pipeline(
        samples, reference, refindex, ref2taxid, taxonomy,
        database, kmer_distribution, template)

    output(results)
}
