include { run_amr } from '../modules/local/amr'
include { run_common } from '../modules/local/common'
include {
    createAbundanceTables;
} from "../modules/local/common"

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process minimap {
    label "wfmetagenomics"
    cpus params.threads
    input:
        tuple val(meta), path(concat_seqs), path(fastcat_stats)
        path reference
        path ref2taxid
        path taxonomy
        val taxonomic_rank
    output:
        tuple(
            val(meta),
            path("*.bam"),
            path("*.bam.bai"),
            emit: bam)
        tuple(
            val(meta),
            path("*assignments.tsv"),
            emit: assignments)
        tuple(
            val(meta),
            path("*lineages.txt"),
            emit: lineage_txt)
        tuple(
            val(meta),
            path("${meta.alias}.json"),
            emit: lineage_json)
    script:
        def sample_id = "${meta.alias}"
        def split = params.split_prefix ? '--split-prefix tmp' : ''
    """
    minimap2 -t "${task.cpus}" ${split} -ax map-ont "${reference}" "${concat_seqs}" \
    | samtools view -h -F 2304 - \
    | workflow-glue format_minimap2 - -o "${sample_id}.minimap2.assignments.tsv" -r "${ref2taxid}" \
    | samtools sort -o "${sample_id}.bam" -
    samtools index "${sample_id}.bam"
    awk -F '\\t' '{print \$3}' "${sample_id}.minimap2.assignments.tsv" > taxids.tmp
    taxonkit \
        --data-dir "${taxonomy}" \
        lineage -R taxids.tmp \
        | workflow-glue aggregate_lineages -p "${sample_id}.minimap2" -r "${taxonomic_rank}"
    file1=`cat *.json`
    echo "{"'"$sample_id"'": "\$file1"}" >> temp
    cp "temp" "${sample_id}.json"
    """
}


process extractMinimap2Reads {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(meta), path("alignment.bam"), path("alignment.bai")
        path ref2taxid
        path taxonomy
    output:
        tuple(
            val(meta),
            path("${meta.alias}.minimap2.extracted.fastq"),
            emit: extracted)
    script:
        def sample_id = "${meta.alias}"
        def policy = params.minimap2exclude ? '--exclude' : ""
    """
    taxonkit \
        --data-dir "${taxonomy}" \
        list -i "${params.minimap2filter}" \
        --indent "" > taxids.tmp
    samtools view -b -F 4 "alignment.bam" > "mapped.bam"
    workflow-glue extract_minimap2_reads \
        "mapped.bam" \
        -r "${ref2taxid}" \
        -o "${sample_id}.minimap2.extracted.fastq" \
        -t taxids.tmp \
        ${policy}
    """
}

// Process to compute the sequencing depth of each reference and their coverages.
process getAlignmentStats {
    label "wfmetagenomics"
    cpus params.threads
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai")
        path ref2taxid
        path taxonomy
    output:
        path "*.tsv.gz", emit: align_stats
    script:
        def sample_name = meta["alias"]
    """
    samtools depth input.bam | bgzip -c > "${sample_name}.depth.tsv.gz" 
    # Reference stats
    samtools coverage input.bam \
    | awk 'BEGIN { OFS="\t" } { if((\$4 != 0) && (\$6 != 0)) {print } }' \
    | sort --parallel=${task.cpus} > "${sample_name}.reference_coverage.tsv"
    # add taxonomy info
    if [ `zcat "${sample_name}.depth.tsv.gz" | head -n 1 | wc -c ` -ne 0 ]
    then
        cut -f1 "${sample_name}.reference_coverage.tsv" | sed '1d'| grep -f - $ref2taxid \
        |  sort --parallel=${task.cpus} | cut -f2 \
        | taxonkit reformat --data-dir $taxonomy -f "{k}\t{K}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -F -I 1 \
        | sed '1 i\\taxid\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' \
        |paste "${sample_name}.reference_coverage.tsv" - \
        | bgzip -c > "${sample_name}.reference.tsv.gz"
        # compress tsv
        bgzip "${sample_name}.reference_coverage.tsv"
    fi

    """
}


process makeReport {
    label "wfmetagenomics"
    input:
        path "read_stats/per-read-stats*.tsv.gz"
        path abundance_table
        path "alignment_stats/*"
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
        def align_stats = params.minimap2_by_reference ? "--align_stats alignment_stats" : ""
        def amr = params.amr as Boolean ? "--amr ${amr}" : ""
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
        --pipeline "minimap" \
        --abundance_threshold "${params.abundance_threshold}"\
        --n_taxa_barplot "${params.n_taxa_barplot}"\
        ${align_stats} \
        ${amr}
    """
}

 
// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wfmetagenomics"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow minimap_pipeline {
    take:
        samples
        reference
        ref2taxid
        taxonomy
        taxonomic_rank
    main:
        outputs = []
        lineages = Channel.empty()
        metadata = samples.map { it[0] }.toList()


        // Run common
        common = run_common(samples)
        software_versions = common.software_versions
        parameters = common.parameters
        samples = common.samples
        // Initial reads QC
        per_read_stats = samples.map {it[2].resolve('per-read-stats.tsv.gz')}.collect()
        // Run Minimap2
  
        mm2 = minimap(
                samples
                | map { [it[0], it[1], it[2] ?: OPTIONAL_FILE ] },
                reference,
                ref2taxid,
                taxonomy,
                taxonomic_rank
            )
        lineages = lineages.mix(mm2.lineage_json)
        // Add some statistics related to the mapping
        if (params.minimap2_by_reference) {
            alignment_reports = getAlignmentStats(mm2.bam, ref2taxid, taxonomy) | collect
        } else {
            alignment_reports = Channel.empty()
        }
        
        abundance_tables = createAbundanceTables(
            lineages.flatMap { meta, lineages_json -> lineages_json }.collect(),
            taxonomic_rank)
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
            alignment_reports.ifEmpty(OPTIONAL_FILE),
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

        if (params.keep_bam) {
            ch_to_publish = ch_to_publish | mix (
            mm2.bam | map { meta, bam, bai -> [bam, "bam"]},
            mm2.bam | map { meta, bam, bai -> [bai, "bam"]},
            )
        }

        // Extract (or exclude) reads belonging (or not) to the chosen taxids.
        if (params.minimap2filter) {
            mm2_filt = extractMinimap2Reads(
                mm2.bam,
                ref2taxid,
                taxonomy
            )
            ch_to_publish = ch_to_publish | mix (
            mm2_filt.extracted | map { meta, fastq -> [fastq, "filtered"]},
            )
        }

        ch_to_publish | output
    emit:
        report.report_html  // just emit something
}


