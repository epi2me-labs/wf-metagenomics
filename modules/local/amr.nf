process abricate{
    label "amr"
    tag "${meta.alias}"
    cpus 1
    memory "8GB"
    input:
        tuple val(meta), path("input_reads.fastq.gz"), path("fastcat_stats/")
        val amr_db
        val amr_minid
        val amr_mincov
    output:
        tuple val(meta), path("${meta.alias}_amr_results.tsv")
    script:
    """
    gunzip -c input_reads.fastq.gz  > input_reads.fastq
    abricate --db $amr_db --minid $amr_minid --mincov $amr_mincov input_reads.fastq > ${meta.alias}_amr_results.tsv
    """
}


process abricate_json{
    label "wfmetagenomics"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path("${meta.alias}_amr_results.tsv")
    output:
        tuple val(meta), path("${meta.alias}_amr.json")
    script:
    """
    workflow-glue abricate_utils \
    --sampleid "${meta.alias}" \
    --input "${meta.alias}_amr_results.tsv" \
    --output "${meta.alias}_amr.json"
    """

}


// Scan step similar to progressive kraken reports 
process progressive_amr{
    label "wfmetagenomics"
    tag "${sample_id}"
    maxForks 1
    cpus 1
    memory "2GB"
    publishDir path: "${params.out_dir}", mode: 'copy', pattern: "${new_state}", saveAs: {name -> "amr"}, overwrite: true
    input:
        path amr_report
        val sample_ids
    output:
        path("${new_state}"), emit: reports
        val(sample_id), emit: sample_id
    script:
        def new_input = amr_report instanceof List ? amr_report.first() : amr_report
        def state = amr_report instanceof List ? amr_report.last() : "NOSTATE"
        sample_id = sample_ids instanceof List ? sample_ids.first() : sample_ids
        new_state = "amr.${task.index}.${sample_id}"
        // n.b where this is used below the files will have been moved, hence new_state
        old_input = "${new_state}/${sample_id}_amr.json"
        """
    if [[ "${task.index}" == "1" ]]; then
        mkdir "${state}"
    fi

    cp -r "${state}" "${new_state}" 
    touch "${old_input}"

    workflow-glue abricate_utils \
        --sampleid $sample_id \
        --input "${new_input}" "${old_input}" \
        --output "${sample_id}_amr.json" \
        --combine
    mv "${sample_id}_amr.json" "${new_state}/${sample_id}_amr.json"
    """
}


workflow run_amr {
    take:
        fastq
        amr_db
        amr_minid
        amr_mincov
    main:
        amr_results = abricate(fastq, amr_db, amr_minid, amr_mincov)
        amr_json = abricate_json(amr_results)
        amr_json.multiMap{ 
            meta, amr ->
            sample_id: meta.alias
            report: amr
        } . set {amr_scan}
        progressive_amr.scan(amr_scan.report, amr_scan.sample_id)
        // amr_all = combine_jsons(amr_json.collect())
    emit:
        reports = progressive_amr.out.reports
}