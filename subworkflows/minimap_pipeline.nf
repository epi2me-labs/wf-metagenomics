import groovy.json.JsonBuilder

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


process combineFilterFastq {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple(
            val(meta.sample_id),
            path("${meta.sample_id}.fastq"),
            emit: filtered)
        tuple(
            val(meta.sample_id),
            path("${meta.sample_id}.stats"),
            emit: stats)
    script:
        def max_length = "${params.max_len}"== null ? "-b ${params.max_len}" : ""
    shell:
    """
    fastcat \
        -a "${params.min_len}" \
        ${max_length} \
        -s "${meta.sample_id}" \
        -r "${meta.sample_id}.stats" \
        -x "${directory}" > "${meta.sample_id}.fastq"
    """
}


process minimap {
    label "wfmetagenomics"
    cpus params.threads
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
            path("${sample_id}.json"),
            emit: lineage_json)
    script:
        def split = params.split_prefix ? '--split-prefix tmp' : ''
    """
    minimap2 -t "${task.cpus}" ${split} -ax map-ont "${reference}" "${reads}" \
    | samtools view -h -F 2304 - \
    | workflow-glue format_minimap2 - -o "${sample_id}.minimap2.assignments.tsv" -r "${ref2taxid}" \
    | samtools sort -o "${sample_id}.bam" -
    samtools index "${sample_id}.bam"
    awk -F '\\t' '{print \$3}' "${sample_id}.minimap2.assignments.tsv" > taxids.tmp
    taxonkit \
        --data-dir "${taxonomy}" \
        lineage -R taxids.tmp \
        | workflow-glue aggregate_lineages -p "${sample_id}.minimap2"
    file1=`cat *.json`
    echo "{"'"$sample_id"'": "\$file1"}" >> temp
    cp "temp" "${sample_id}.json"


    """
}


process extractMinimap2Reads {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path("alignment.bam"), path("alignment.bai")
        path ref2taxid
        path taxonomy
    output:
        tuple(
            val(sample_id),
            path("${sample_id}.minimap2.extracted.fastq"),
            emit: extracted)
    script:
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
        path "lineages/*"
        path "versions/*"
        path "params.json"
        path template
    output:
        path "wf-metagenomics-*.html", emit: report_html
    script:
        report_name = "wf-metagenomics-report.html"
    """
    workflow-glue report \
        "${report_name}" \
        --versions versions \
        --params params.json \
        --summaries ${stats} \
        --lineages lineages \
        --vistempl "${template}" \
        --pipeline "minimap"
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
workflow minimap_pipeline {
    take:
        samples
        reference
        refindex
        ref2taxid
        taxonomy
        template
    main:
        outputs = []
        lineages = Channel.empty()
        taxonomy = unpackTaxonomy(taxonomy)
    
        // Initial reads QC
        reads = combineFilterFastq(samples)

        // Run Minimap2
  
        mm2 = minimap(
                reads.filtered,
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

        output(report.report_html.mix(
            software_versions, workflow_params))
    emit:
        report.report_html  // just emit something
}


