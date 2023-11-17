import nextflow.util.BlankSeparatedList

include { run_amr } from '../modules/local/amr'
include {
    run_common;
    createAbundanceTables;
    output;
} from "../modules/local/common"


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

// Filter reads, calculate some stats, and run kraken2 classification
process run_kraken2 {
    label 'wfmetagenomics'
    publishDir "${params.out_dir}/kraken2", mode: 'copy', pattern: "*kraken2*"
    cpus 1
    input:
        tuple val(meta), path(sample_fastq), path(fastq_stats)
        path kraken_db
    output:
        tuple val(meta), path("${meta.alias}.kraken2.report.txt"), path("${meta.alias}.kraken2.assignments.tsv"), emit: kraken2_reports
    script:
        def sample_id = "${meta.alias}"
        def memory_mapping = params.kraken2_memory_mapping ? '--memory-mapping' : ''
    """
    kraken2 --db ${kraken_db} ${sample_fastq}\
    --threads $params.threads \
    --report "${sample_id}.kraken2.report.txt"\
    --confidence ${params.kraken2_confidence} ${memory_mapping} > "${sample_id}.kraken2.assignments.tsv"
    """
}


process run_bracken {
    label "wfmetagenomics"
    publishDir "${params.out_dir}/bracken", mode: 'copy', pattern: "*bracken*"
    cpus 2
    input:
        tuple val(meta), path("kraken2.report"), path("kraken2.assignments.tsv")
        path(database)
        path(taxonomy)
        path "bracken_length.txt"
        val(taxonomic_rank)
    output:
        tuple val(meta), path("${meta.alias}.kraken2_bracken.report"), emit: bracken_reports
        tuple val(meta), path("${meta.alias}.json"), emit: bracken_json
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
    
    # add sample to the json file    
    file1=`cat *.json`
    echo "{"'"$sample_id"'": "\$file1"}" >> "bracken.json"
    mv "bracken.json" "${sample_id}.json"
    """
}

// Concatenate kraken reports per read
process output_kraken2_read_assignments {
    label "wfmetagenomics"
    input:
        tuple val(meta), path("${meta.alias}.kraken2.assignments.tsv")
        path taxonomy
    output:
        tuple(
            val(meta),
            path("*_lineages.kraken2.assignments.tsv"), 
            emit: kraken2_reads_assignments
            )
    script:
        def sample_id = "${meta.alias}"
    """
    # Run taxonkit to give users a more informative table
    taxonkit reformat  -I 3  --data-dir "${taxonomy}" -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}" -F "${sample_id}.kraken2.assignments.tsv" > "${sample_id}_lineages.kraken2.assignments.tsv"
    """

}

process makeReport {
    label "wfmetagenomics"
    input:
        path "read_stats/per-read-stats*.tsv.gz"
        path abundance_table
        path "lineages/*"
        path "versions/*"
        path "params.json"
        val taxonomic_rank
        path amr
    output:
        path "*.html", emit: report_html
    script:
        String workflow_name = workflow.manifest.name.replace("epi2me-labs/","")
        String report_name = "${workflow_name}-report.html"
        def stats_args = params.wf.stats ? "--read_stats read_stats/*" : ""
        amr = params.amr as Boolean ? "--amr ${amr}" : ""
    """
    workflow-glue report \
        "${report_name}" \
        --workflow_name ${workflow_name} \
        --versions versions \
        --params params.json \
        ${stats_args} \
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
    main:
        OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

        // Run common
        common = run_common(samples)
        software_versions = common.software_versions
        parameters = common.parameters
        samples = common.samples
        // Initial reads QC
        per_read_stats = samples.map {it[2].resolve('per-read-stats.tsv.gz')}.collect()
        | ifEmpty ( OPTIONAL_FILE )
        // Run Kraken2
        kraken2_reports = run_kraken2(samples, database)
        
        // Run bracken
        bracken_reports = run_bracken(kraken2_reports, database, taxonomy, bracken_length, taxonomic_rank)
        lineages = bracken_reports.bracken_json

        // Abundance tabeles
        abundance_tables = createAbundanceTables(
            lineages.flatMap { meta, lineages_json -> lineages_json }.collect(),
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
            per_read_stats,
            abundance_tables.abundance_tsv,
            lineages.flatMap { meta, lineages_json -> lineages_json }.collect(),
            software_versions,
            parameters,
            taxonomic_rank,
            amr_reports.ifEmpty(OPTIONAL_FILE)
        )

        ch_to_publish = Channel.empty()
        | mix(
            software_versions,
            parameters,
            report.report_html,
            abundance_tables.abundance_tsv,
        )
        | map { [it, null] }

        // output kraken read assignments + taxonomy info
        if (params.include_kraken2_assignments) {
            kraken2_assignments = output_kraken2_read_assignments(
                kraken2_reports.map{ 
                    id, report, assignments -> tuple(id, assignments) 
                }.groupTuple(), taxonomy)
            ch_to_publish = ch_to_publish | mix (
            kraken2_assignments.kraken2_reads_assignments | map {
                 id, kraken_reads_classification -> [kraken_reads_classification, "kraken_reads_assignments"]},
            )
        }

         ch_to_publish | output

    emit:
        report.report_html  // just emit something
}

