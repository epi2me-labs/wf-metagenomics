import groovy.json.JsonBuilder

include { configure_igv } from '../lib/common'
include { run_amr } from '../modules/local/amr'
include { filter_references } from '../modules/local/igv_related'
include {
    run_common;
    createAbundanceTables;
    publish;
    publishReads;
} from "../modules/local/common"

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process minimap {
    label "wfmetagenomics"
    tag "${meta.alias}"
    cpus params.threads
    publishDir (
        "${params.out_dir}/bams", mode: 'copy',
        pattern: "${meta.alias}.reference.bam*",
        enabled: (params.keep_bam || params.igv))
    publishDir (
        "${params.out_dir}/bams", mode: 'copy',
        pattern: "${meta.alias}.bamstats_results",
        enabled: (params.keep_bam || params.igv))
    publishDir (
        "${params.out_dir}/reads_assignments_unfiltered", mode: 'copy',
        pattern: "${meta.alias}.minimap2.assignments.tsv",
        enabled: params.include_read_assignments)
    // due to the wf fail at the samtools step
    memory {12.GB * task.attempt}
    maxRetries 1
    errorStrategy = 'retry'
    input:
        tuple val(meta),
            path(concat_seqs),
            path(stats)
        path reference
        path ref2taxid
        val common_minimap2_opts
    output:
        tuple val(meta),
            path("${meta.alias}.reference.bam"),
            path("${meta.alias}.reference.bam.bai"),
            path("${meta.alias}.bamstats_results"),
            env(n_unmapped),
            emit: bam
        tuple val(meta),
            path("${meta.alias}.modified.minimap2.assignments.tsv"),
            emit: assignments
        tuple val(meta),
            path("${meta.alias}.minimap2.assignments.tsv"),
            emit: raw_assignments, optional:true
    script:
        def common_minimap2_opts = (reference.size() > 4e9 ) ? common_minimap2_opts + ['--split-prefix tmp'] : common_minimap2_opts
        common_minimap2_opts = common_minimap2_opts.join(" ")
        def bamstats_threads = Math.max(1, task.cpus - 1)
    // min_percent_identity and min_ref_coverage can be used within format_minimap2 or after the BAM is generated to not modify the raw BAM from the alignment.
    // Filter from ${meta.alias}.minimap2.assignments.tsv
    """
    minimap2 -t $task.cpus ${common_minimap2_opts} $reference $concat_seqs \
    | samtools view -h -F 2304 - \
    | workflow-glue format_minimap2 - -o "${meta.alias}.minimap2.assignments.tsv" -r "$ref2taxid" \
    | samtools sort -@ ${task.cpus - 1} --write-index -o "${meta.alias}.reference.bam##idx##${meta.alias}.reference.bam.bai" -

    # run bamstats
    mkdir "${meta.alias}.bamstats_results"
    bamstats "${meta.alias}.reference.bam" -s $meta.alias -u \
        -f ${meta.alias}.bamstats_results/${meta.alias}.bamstats.flagstat.tsv -t $bamstats_threads \
        --histograms histograms \
    | bgzip > ${meta.alias}.bamstats_results/bamstats.readstats.tsv.gz
    rm -r histograms/

    # Get readIDs that do not satisfy the minimum percentage of identity or ref_cov filters.
    csvtk filter2 -t "${meta.alias}.bamstats_results/bamstats.readstats.tsv.gz" -f '\$ref_coverage < ${params.min_ref_coverage} || \$iden < ${params.min_percent_identity}' \
        | csvtk cut -t -f name - \
        | sed '1d' - > extra_unclassified_readsIDS.txt
    # if there are reads to be reclassified
    if [ -s extra_unclassified_readsIDS.txt ]; then
        # Reclassified minimap2.assignments.tsv and move to unclassified those reads which haven't passed the filter
        grep -w -f extra_unclassified_readsIDS.txt -F "${meta.alias}.minimap2.assignments.tsv" \
            | csvtk replace -t -H -f 3 -p '(.+)' -r 0 - > modified.minimap2.assignments.txt
        # Avoid grep error when there are no hits and fails with set -eo pipefail
        { grep -w -v -f extra_unclassified_readsIDS.txt -F "${meta.alias}.minimap2.assignments.tsv" || true; }  \
            | cat - modified.minimap2.assignments.txt > "${meta.alias}.modified.minimap2.assignments.tsv"
    else
        cp "${meta.alias}.minimap2.assignments.tsv" "${meta.alias}.modified.minimap2.assignments.tsv"
    fi

    # get unmapped reads to add it to meta (unmapped == unclassified)
    n_unmapped=\$(csvtk grep -t -f ref -p '*' "${meta.alias}.bamstats_results/${meta.alias}.bamstats.flagstat.tsv" \
        | csvtk cut -t -f total - \
        | sed "1d" -
    )
    """
}


process minimapTaxonomy {
    label "wfmetagenomics"
    tag "${meta.alias}"
    cpus 1
    publishDir "${params.out_dir}/reads_assignments", mode: 'copy', pattern: "*_lineages.minimap2.assignments.tsv", enabled: params.include_read_assignments
    memory 4.GB
    input:
        tuple val(meta), path(assignments)
        path taxonomy
        val taxonomic_rank
    output:
        tuple val(meta),
            path("${meta.alias}.json"),
            emit: lineage_json
        tuple val(meta),
            path("${meta.alias}.unclassified.txt"),
            emit: unclassified_ids
        tuple val(meta),
            path("${meta.alias}_lineages.minimap2.assignments.tsv"),
            emit: reads_assignments, optional:true
    script:
    String taxid_col = 3
    """
    # add taxonomy to the tables if required by the user,
    # as this file would contain 1 entry per read
    # Output assignments after filters
    if ${params.include_read_assignments} ;
    then
        taxonkit reformat \
            --taxid-field $taxid_col \
            --data-dir "${taxonomy}" \
            --format "{k}|{p}|{c}|{o}|{f}|{g}|{s}" \
            --fill-miss-rank "${assignments}" > "${meta.alias}_lineages.minimap2.assignments.tsv"
    fi

    awk '{print \$3}' "${assignments}" > taxids.tmp
    # Add taxonomy
    taxonkit \
        --data-dir "${taxonomy}" \
        lineage --show-lineage-ranks taxids.tmp \
        | workflow-glue aggregate_lineages -p "${meta.alias}.minimap2" -r "${taxonomic_rank}"
    file1=\$(find -name '*.json' -exec cat {} +)
    echo "{"'"${meta.alias}"'": \$file1}" >> temp
    cp "temp" "${meta.alias}.json"

    # Recover unclassified IDs
    csvtk filter2 --no-header-row --tabs -f '\$1=="U"' "${assignments}" \
        | cut -f2 > "${meta.alias}.unclassified.txt"
    """
}


process extractMinimap2Reads {
    label "wfmetagenomics"
    publishDir "${params.out_dir}/extracted", mode: 'copy', pattern: "*.minimap2.extracted.fastq"
    tag "${meta.alias}"
    cpus 1
    memory "7 GB" //depends on the size of the BAM file.
    input:
        tuple val(meta),
            path("alignment.bam"),
            path("alignment.bai"),
            path("bamstats"),
            val(n_unmapped)
        path ref2taxid
        path taxonomy
    output:
        tuple val(meta),
            path("${meta.alias}.minimap2.extracted.fastq"),
            emit: extracted
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


/* Process to compute the sequencing depth of each reference and their coverages.
Run python script to parse the data and output the alignment table
*/
process getAlignmentStats {
    label "wfmetagenomics"
    tag "${meta.alias}"
    publishDir "${params.out_dir}/alignment_tables", mode: 'copy', pattern: "*-alignment-stats.tsv"
    cpus Math.max(params.threads, 2)
    //depends on number of references and their lengths. There are also custom databases of varying sizes.
    memory "7 GB"
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai"), path("bamstats")
        path ref2taxid
        path taxonomy
    output:
        path "${meta.alias}-alignment-stats*", emit: align_stats
    script:
        def sample_name = meta["alias"]
    // TODO: remove samtools coverage and use bamstats results
    """
    samtools depth input.bam | bgzip -c > "${sample_name}.depth.tsv.gz"
    # Reference stats
    samtools coverage input.bam \
    | awk 'BEGIN { OFS="\t" } { if((\$4 != 0) && (\$6 != 0)) {print } }' \
    | sort --parallel=${task.cpus - 1} > "${sample_name}.reference_coverage.tsv"
    # add taxonomy info
    if [ `zcat "${sample_name}.depth.tsv.gz" | head -n 1 | wc -c ` -ne 0 ]
    then
        cut -f1 "${sample_name}.reference_coverage.tsv" | sed '1d'| grep -w -f - $ref2taxid \
        |  sort --parallel=${task.cpus - 1} | cut -f2 \
        | taxonkit reformat --data-dir $taxonomy -f "{k}\t{K}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}" -F -I 1 \
        | sed '1 i\\taxid\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' \
        |paste "${sample_name}.reference_coverage.tsv" - \
        | bgzip -c > "${sample_name}.reference.tsv.gz"
        # compress tsv
        bgzip "${sample_name}.reference_coverage.tsv"

        # Run python script to process a useful table
        workflow-glue alignment_stats \
        --output "${meta.alias}-alignment-stats.tsv" \
        --output_heatmap "${meta.alias}-alignment-stats-heatmap.tsv" \
        --coverage "${sample_name}.reference.tsv.gz" \
        --depth "${sample_name}.depth.tsv.gz"
    fi
    """
}


process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "${report_name}"
    cpus 1
    memory "4 GB"
    input:
        val wf_version
        val metadata
        path(stats, stageAs:"stats/stats_*")
        path abundance_table
        path "alignment_tables/*"
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
        String align_stats = params.minimap2_by_reference ? "--align_stats alignment_tables" : ""
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
        --pipeline "minimap2" \
        --abundance_threshold "${params.abundance_threshold}"\
        --n_taxa_barplot "${params.n_taxa_barplot}"\
        ${align_stats} \
        ${amr}
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
        common_minimap2_opts
        keep_bam
        output_igv
    main:
        lineages = Channel.empty()
        // Run common
        common = run_common(
            samples,
            common_minimap2_opts)
        software_versions = common.software_versions
        parameters = common.parameters
        samples = common.samples

        // Run Minimap2
        mm2 = minimap(
                samples,
                reference,
                ref2taxid,
                common_minimap2_opts
        )
        // Add taxonomy
        mm2_taxonomy = minimapTaxonomy(mm2.assignments, taxonomy, taxonomic_rank)
        // add unclassified to meta to use it to filter samples with all unclassified in igv
        samples_classification = mm2.bam.map { meta, bam, bai, stats, unmapped->
            [meta + [n_unclassified: unmapped as Integer], bam, bai, stats]
        }


        // Output unclassified
        if (params.output_unclassified) {
            unclassified_to_extract = samples.join(
                    mm2_taxonomy.unclassified_ids
                ).map { meta, seqs, stats, unclassified_ids ->
                        [meta, seqs, unclassified_ids]
                }
            publishReads(unclassified_to_extract, "unclassified")
        }
        // Use initial reads stats (after fastcat) QC, but update meta
        for_report = samples
        | map{
            meta, path, stats -> [meta.alias, stats]
        }
        | combine (
            samples_classification
            | map {
                    meta, bam, bai, stats -> [meta.alias, meta]
                },
            by: 0
        ) | multiMap{ alias, stats, meta ->
            meta: meta
            stats: stats
        }
        metadata = for_report.meta.collect()
        // create a file list of the stats, and signal if its empty or not
        stats = for_report.stats.collect()
        // take samples with classified sequences
        bam_classified = samples_classification.filter { meta, bam, bai, stats ->
            // a sample is unclassified if all reads are unclassified
            meta.n_seqs > meta.n_unclassified
        }

        lineages = lineages.mix(mm2_taxonomy.lineage_json)
        // Add some statistics related to the mapping
        if (params.minimap2_by_reference) {
            alignment_reports = getAlignmentStats(bam_classified, ref2taxid, taxonomy) | collect
        } else {
            alignment_reports = Channel.empty()
        }
        abundance_tables = createAbundanceTables(
            lineages.flatMap { meta, lineages_json -> lineages_json }.collect(),
            taxonomic_rank, params.classifier)

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
            alignment_reports.ifEmpty(OPTIONAL_FILE),
            lineages.flatMap { meta, lineages_json -> lineages_json }.collect(),
            software_versions,
            parameters,
            taxonomic_rank,
            amr_reports.ifEmpty(OPTIONAL_FILE)
        )

        if (output_igv) {
            // filter references
            bamstats_flagstat = bam_classified
            | map { meta, bam, bai, stats -> file(stats.resolve("*.bamstats.flagstat.tsv")) }
            | collect

            filtered_refs = filter_references(
                reference,
                bamstats_flagstat
            )
            // create IGV config file: use absolute paths in the igv.JSON
            // write files in file-names.txt

            igv_files = filtered_refs
                | map { list -> list.collect { "igv_reference/${it.Name}" }}
                | concat (
                    bam_classified
                    | flatMap {
                        meta, bam, bai, stats -> [
                            "bams/${meta.alias}.reference.bam",
                            "bams/${meta.alias}.reference.bam.bai",
                        ]
                    }
                    | collect
                )
                | flatten
                | collectFile(name: "file-names.txt", newLine: true, sort: false)
            boolean keep_track_order = false
            igv_conf = configure_igv(
                igv_files,
                Channel.of(null), // igv locus
                [displayMode: "SQUISHED", colorBy: "strand"], // bam extra opts
                Channel.of(null), // vcf extra opts
                keep_track_order
                )
        }


        // Extract (or exclude) reads belonging (or not) to the chosen taxids.
        if (params.minimap2filter) {
            mm2_filt = extractMinimap2Reads(
                mm2.bam,
                ref2taxid,
                taxonomy
            )
        }

    emit:
        report.report_html  // just emit something
}


