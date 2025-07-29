process abricate{
    label "amr"
    tag "${meta.alias}"
    cpus 2
    memory {7.GB * task.attempt}
    maxRetries 3
    errorStrategy 'retry'
    input:
        tuple val(meta), path("input_reads.fastq.gz"), path("stats/")
        val amr_db
        val amr_minid
        val amr_mincov
    output:
        tuple val(meta), path("${meta.alias}_amr_results.tsv")
    script:
        String fastq_name = "${meta.alias}.fastq"
    """
    # run a patched (this is done in the container) version of abricate, which uses `seqkit fq2fa` in place of `any2fasta`
    abricate --db $amr_db --minid $amr_minid --mincov $amr_mincov input_reads.fastq.gz --threads $task.cpus > "${meta.alias}_amr_results.tsv"
    """
}


process abricate_json{
    label "wfmetagenomics"
    publishDir path: "${params.out_dir}/amr", mode: 'copy', pattern: "*.json", saveAs: { filename -> "${meta.alias}.amr.json" }
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path("${meta.alias}_amr_results.tsv")
    output:
        tuple val(meta), path("${meta.alias}_amr_per_chunk.json")
    script:
    """
    workflow-glue abricate_utils \
    --sampleid "${meta.alias}" \
    --input "${meta.alias}_amr_results.tsv" \
    --output "${meta.alias}_amr_per_chunk.json"
    """

}


workflow run_amr {
    take:
        input
        amr_db
        amr_minid
        amr_mincov
    main:
        amr_results = abricate(input, amr_db, amr_minid, amr_mincov)
        amr_json = abricate_json(amr_results)
        amr_all = amr_json.flatMap { meta, amr_json -> amr_json }.collect()
    emit:
        reports = amr_all
}