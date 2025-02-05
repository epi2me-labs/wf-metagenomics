import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList

include { run_amr } from '../modules/local/amr'
include {
    run_common;
    createAbundanceTables;
    publish;
} from "../modules/local/common"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// Filter reads, calculate some stats, and run kraken2 classification
process run_kraken2 {
    label 'wfmetagenomics'
    tag "${meta.alias}"
    publishDir "${params.out_dir}/kraken2", mode: 'copy', pattern: "*kraken2.report.txt*"
    cpus params.threads
    // Set the memory required to the size of the database + 4GB overhead.
    memory {
        if (params.kraken2_memory_mapping) {
            "${hash_size.intdiv(4) + 1e9} B " // use less memory if memory_mapping is used
        } else {
            "${hash_size + 4e9} B "
        }
    }
    errorStrategy {
        task.exitStatus == 137 ? log.error("Error 137 may indicate the process ran out of memory.\nIf you are using Docker you should check the amount of RAM allocated to your Docker server.") : ''
        log.error("Consider to use --kraken2_memory_mapping to reduce the use of RAM memory.")
    }
    input:
        tuple val(meta), path(concat_seqs), path(fastq_stats)
        path kraken_db
        val hash_size
    output:
        tuple val(meta), path("${meta.alias}.kraken2.report.txt"), path("${meta.alias}.kraken2.assignments.tsv"), emit: kraken2_reports
    script:
        def sample_id = "${meta.alias}"
        def memory_mapping = params.kraken2_memory_mapping ? '--memory-mapping' : ''
    """
    kraken2 --db ${kraken_db} ${concat_seqs} \
        --threads $task.cpus \
        --report "${sample_id}.kraken2.report.txt" \
        --confidence ${params.kraken2_confidence} ${memory_mapping} > "${sample_id}.kraken2.assignments.tsv"
    """
}


process run_bracken {
    label "wfmetagenomics"
    tag "${meta.alias}"
    publishDir "${params.out_dir}/bracken", mode: 'copy', pattern: "*bracken*"
    cpus Math.max(params.threads - 2, 2)
    memory {8.GB * task.attempt - 1.GB}
    maxRetries 1
    errorStrategy = 'retry'
    input:
        tuple val(meta), path("kraken2.report"), path("kraken2.assignments.tsv")
        path(database)
        path(taxonomy)
        path "bracken_length.txt"
        val(taxonomic_rank)
    output:
        tuple val(meta), path("${meta.alias}.kraken2_bracken.report"), emit: bracken_reports
        tuple val(meta), path("${meta.alias}.json"), env(n_unclassified), emit: bracken_json
    script:
        def sample_id = "${meta.alias}"
        def awktab="awk -F '\t' -v OFS='\t'"
    """
    # run bracken on the latest kreports, is this writing some outputs
    # alongside the inputs? seems at least {}.kreport_bracken_species.txt
    # is written alongside the input
    BRACKEN_LENGTH=\$(cat bracken_length.txt)

    workflow-glue run_bracken \
        "${database}" \
        kraken2.report \
        \$BRACKEN_LENGTH \
        "${taxonomic_rank}" \
        "${params.bracken_threshold}" \
        "${sample_id}.kraken2_bracken.report"

    # do some stuff...
    ${awktab} '{ print \$2,\$6 }' "${sample_id}.kraken2_bracken.report" \
    | ${awktab} 'NR!=1 {print}' \
    | tee taxacounts.txt \
    | ${awktab} '{ print \$1 }' > taxa.txt
    taxonkit lineage \
        -j ${task.cpus} \
        --data-dir $taxonomy \
        -R taxa.txt  > lineages.txt
    workflow-glue aggregate_lineages_bracken \
        -i "lineages.txt" -b "taxacounts.txt" \
        -u kraken2.report \
        -p "${sample_id}.kraken2" \
        -r "${taxonomic_rank}"
    n_unclassified=\$(cut -f1 kraken2.assignments.tsv | { grep -c '^U' - || test \$? = 1;} )
    # add sample to the json file    
    file1=\$(find -name '*.json' -exec cat {} +)
    echo "{"'"$sample_id"'": \$file1}" >> "bracken.json"
    mv "bracken.json" "${sample_id}.json"
    """
}

// Concatenate kraken reports per read
process output_kraken2_read_assignments {
    label "wfmetagenomics"
    publishDir "${params.out_dir}/reads_assignments", mode: 'copy', pattern: "*_lineages.kraken2.assignments.tsv"
    tag "${meta.alias}"
    cpus 2
    memory "4 GB"
    input:
        tuple val(meta), path("${meta.alias}.kraken2.assignments.tsv")
        path taxonomy
    output:
        tuple(
            val(meta),
            path("*_lineages.kraken2.assignments.tsv")
            )
    script:
        def sample_id = "${meta.alias}"
    """
    # Run taxonkit to give users a more informative table
    taxonkit reformat \
        --taxid-field 3 \
        --data-dir "${taxonomy}" \
        --format "{k}|{p}|{c}|{o}|{f}|{g}|{s}" \
        --fill-miss-rank "${sample_id}.kraken2.assignments.tsv" > "${sample_id}_lineages.kraken2.assignments.tsv"
    """

}

process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "${report_name}"
    cpus 1
    // Report generation will generally use less memory than 4GB, but higher complexity data will use more.
    memory {4.GB * task.attempt}
    maxRetries 3
    errorStrategy = 'retry'
    input:
        val wf_version
        val metadata
        path(stats, stageAs:"stats/stats_*")
        path abundance_table
        path "lineages/*"
        path "versions/*"
        path "params.json"
        val taxonomic_rank
        path "amr/*"
    output:
        path "${report_name}", emit: report_html
    script:
        String workflow_name = workflow.manifest.name.replace("epi2me-labs/","")
        String metadata = new JsonBuilder(metadata).toPrettyString()
        report_name = "${workflow_name}-report.html"
        String amr = params.amr as Boolean ? "--amr amr" : ""
    """
    echo '${metadata}' > metadata.json
    workflow-glue report \
        "${report_name}" \
        --workflow_name ${workflow_name} \
        --versions versions \
        --params params.json \
        --wf_version $wf_version \
        --metadata metadata.json \
        --read_stats $stats \
        --lineages lineages \
        --abundance_table "${abundance_table}" \
        --taxonomic_rank "${taxonomic_rank}" \
        --pipeline "kraken2" \
        --abundance_threshold "${params.abundance_threshold}"\
        --n_taxa_barplot "${params.n_taxa_barplot}"\
        ${amr}
    """
}

workflow kraken_pipeline {
    take:
        samples
        taxonomy
        database
        bracken_length
        taxonomic_rank
        common_minimap2_opts

    main:
        OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
        // Run common
        common = run_common(samples, common_minimap2_opts)
        software_versions = common.software_versions
        parameters = common.parameters
        samples = common.samples

        // Run Kraken2
        // Find out size of the db. Cannot be done within the process
        database_main_file_size = database.resolve('hash.k2d').size()
        kraken2_reports = run_kraken2(samples, database, database_main_file_size)

        // Run bracken
        bracken_reports = run_bracken(kraken2_reports, database, taxonomy, bracken_length, taxonomic_rank)
        lineages = bracken_reports.bracken_json
        // Update meta with unclassified
        samples_classification = lineages.map { meta, lineages_json, n_unclassified->
            [meta + [n_unclassified: n_unclassified as Integer], lineages_json]
        }

        // // Use initial reads stats (after fastcat) QC, but update meta
        for_report = samples
        | map{
            meta, path, stats -> [meta.alias, stats]
        }
        | combine (
            samples_classification
            | map {
                    meta, lineages_json -> [meta.alias, meta]
                },
            by: 0
        ) | multiMap{ alias, stats, meta ->
            meta: meta
            stats: stats
        }
        metadata = for_report.meta.collect()
        // create a file list of the stats, and signal if its empty or not
        stats = for_report.stats.collect()

        // Abundance table
        abundance_tables = createAbundanceTables(
            lineages.flatMap { meta, lineages_json, n_unclassified -> lineages_json }.collect(),
            taxonomic_rank, 'kraken2')


        // Process AMR
        if (params.amr) {
            run_amr = run_amr(
                samples,
                "${params.amr_db}",
                "${params.amr_minid}",
                "${params.amr_mincov}"
            )
            amr_reports = run_amr.reports
        } else {
            amr_reports = Channel.empty()
        }

        // Reporting
        report = makeReport(
            workflow.manifest.version,
            metadata,
            stats,
            abundance_tables.abundance_tsv,
            lineages.flatMap { meta, lineages_json, n_unclassified -> lineages_json }.collect(),
            software_versions,
            parameters,
            taxonomic_rank,
            amr_reports.ifEmpty(OPTIONAL_FILE)
        )


        // output kraken read assignments + taxonomy info
        if (params.include_read_assignments) {
            kraken2_assignments = output_kraken2_read_assignments(
                kraken2_reports.map{
                    id, report, assignments -> tuple(id, assignments)
                }.groupTuple(), taxonomy)
        }


    emit:
        report.report_html  // just emit something
}

