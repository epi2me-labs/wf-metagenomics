process abricate{
    label "amr"
    tag "${meta.alias}"
    cpus 1
    memory "7GB"
    input:
        tuple val(meta), path(concat_seqs), path("stats/")
        val amr_db
        val amr_minid
        val amr_mincov
    output:
        tuple val(meta), path("${meta.alias}_amr_results.tsv")
    script:
        String fastq_name = "${meta.alias}.fastq"
    """
    # run abricate
    gunzip -c ${concat_seqs}  > input_reads.fastq
    abricate --db $amr_db --minid $amr_minid --mincov $amr_mincov input_reads.fastq > ${meta.alias}_amr_results.tsv
    """
}


process abricate_json{
    label "wfmetagenomics"
    publishDir path: "${params.out_dir}/amr", mode: 'copy', pattern: "*.json", saveAs: { filename -> "${meta.alias}.amr.json" }, enabled: !params.real_time
    cpus 1
    memory "2GB"
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


// Scan step similar to progressive kraken reports
process progressive_amr{
    label "wfmetagenomics"
    tag "${current_sample_id}"
    maxForks 1
    cpus 1
    memory "2GB"
    publishDir path: "${params.out_dir}", mode: 'copy', pattern: "${new_state}", saveAs: {name -> "amr"}, overwrite: true
    input:
        path report_and_prev_states
        val sample_ids
    output:
        // new state is a folder that contains a JSON file per sample.
        // this is taken as input in the scan process after the new JSON file to be processed (new_input_sample).
        // so if a sample is seen >1: need to combine the new JSON (new_input_sample)
        // with the existing one (old_sample_input)
        // at the same time, next iter is also the last element of the input data: previous_state
        path(new_state), emit: reports
        val(current_sample_id), emit: sample_id
    script:
        // new_input_sample is the first element of the list of inputs, is the "new" element
        // (a JSON file coming from the previous process abricate_json)
        // then each element is a directory, last of them is previous state.
        // previous state is a directory that includes the latests JSON file converted from TSV to JSON
        // the name of this new_input_sample is "${current_sample_id}_amr.json"
        def new_input_sample = report_and_prev_states instanceof List ? report_and_prev_states.first() : report_and_prev_states
        def previous_state = report_and_prev_states instanceof List ? report_and_prev_states.last() : "NOSTATE"
        current_sample_id = sample_ids instanceof List ? sample_ids.first() : sample_ids
        new_state = "amr.${task.index}.${current_sample_id}"
        // old_sample_input doesn't exist first time a sample is processed.
        // in the next iterations (chunks of the same sample) it exists because has been created the first time
        // n.b where this is used below the files will have been moved, hence new_state
        old_sample_input = "${new_state}/${current_sample_id}_amr.json" // this is copied from previous state into new_state
        """
        if [[ "${task.index}" == "1" ]]; then
            mkdir "${previous_state}"
        fi
        cp -r "${previous_state}" "${new_state}"

        # if there are previous outputs of the samples
        # combine them
        # otherwise there is no point in combining anything, new_sample_input is the only input
        # we rename it to the standard final JSON file
        if [[ -f "${old_sample_input}" ]]
            then
            # add info of file old_sample_input to the corresponding file in the dir new_input
            workflow-glue abricate_utils \
                --sampleid "${current_sample_id}" \
                --input "${new_input_sample}" "${old_sample_input}" \
                --output "${current_sample_id}_combined_amr.json" \
                --combine
            mv "${current_sample_id}_combined_amr.json" "${new_state}/${current_sample_id}_amr.json"
        else
            mv "${new_input_sample}" "${new_state}/${current_sample_id}_amr.json"
        fi
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
        amr_json.multiMap{ 
            meta, amr ->
            sample_id: meta.alias
            report: amr
        } . set {amr_scan}
        if (params.real_time) {
            progressive_amr.scan(amr_scan.report, amr_scan.sample_id)
        } else {
            amr_all = amr_json.flatMap { meta, amr_json -> amr_json }.collect()
        }
    emit:
        reports = (params.real_time) ? progressive_amr.out.reports : amr_all
}