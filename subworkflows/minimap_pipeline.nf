include { run_amr } from '../modules/local/amr'
include { run_common } from '../modules/local/common'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

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
    elif [[ "${taxonomy}" == *.zip ]]
    then
        mkdir taxonomy_dir
        unzip "${taxonomy}" -d taxonomy_dir
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


process check_reference_ref2taxid {
    label "wfmetagenomics"
    input:
        path reference
        path ref2taxid
    output:
        val true
    """
    # just check if the reference is a fasta file. It can be a MMI index file.
    if [[ "${reference}" != *.mmi ]]
    then
        if [ "\$(wc -l < "${ref2taxid}")" -eq "\$(grep '>' < "${reference}" | wc -l)" ]
        then
            echo 'Match!'
        else
            echo "Error: The number of elements of the "${ref2taxid}" doesn't match the number of elements in the "${reference}"."
            echo "Please provide the fitting "${ref2taxid}" for the "${reference}"."
            echo "Exiting".
            exit 1
        fi
    fi
    """
}


process minimap {
    label "wfmetagenomics"
    cpus params.threads
    input:
        tuple val(meta), path(concat_seqs), path(fastcat_stats)
        path reference
        path refindex
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


process makeReport {
    label "wfmetagenomics"
    input:
        path per_read_stats
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
        def stats_args = (per_read_stats.name == OPTIONAL_FILE.name) ? "" : "--stats $per_read_stats"
        amr = params.amr as Boolean ? "--amr ${amr}" : ""
    """   
    workflow-glue report \
        "${report_name}" \
        --workflow_name ${workflow_name} \
        --versions versions \
        --params params.json \
        ${stats_args} \
        --lineages lineages \
        --taxonomic_rank "${taxonomic_rank}" \
        --pipeline "minimap" \
        --abundance_threshold "${params.abundance_threshold}"\
        --n_taxa_barplot "${params.n_taxa_barplot}"\
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
        refindex
        ref2taxid
        taxonomy
        taxonomic_rank
    main:
        outputs = []
        lineages = Channel.empty()
        taxonomy = unpackTaxonomy(taxonomy)
        // Initial reads QC
        per_read_stats = samples.map {
            it[2] ? it[2].resolve('per-read-stats.tsv') : null
        }
        | collectFile ( keepHeader: true )
        | ifEmpty ( OPTIONAL_FILE )
        metadata = samples.map { it[0] }.toList()

        // check that the number of elements of the reference (if it is a fasta)
        // matches the number of elements of the ref2taxid.
        check_reference_ref2taxid(reference, ref2taxid)

        // Run Minimap2
  
        mm2 = minimap(
                samples
                | map { [it[0], it[1], it[2] ?: OPTIONAL_FILE ] },
                reference,
                refindex,
                ref2taxid,
                taxonomy,
                taxonomic_rank
            )
        lineages = lineages.mix(mm2.lineage_json)

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
        common = run_common()
        software_versions = common.software_versions
        parameters = common.parameters
        report = makeReport(
            per_read_stats,
            lineages.flatMap { it -> [ it[1] ] }.collect(),
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


