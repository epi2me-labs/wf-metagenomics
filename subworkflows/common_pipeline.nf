process exclude_host_reads {
    label "wfmetagenomics"
    publishDir "${params.out_dir}/host_bam", mode: 'copy', pattern: "*.host.bam*"
    publishDir "${params.out_dir}/no_host_bam", mode: 'copy', pattern: "*.unmapped.bam*"
    tag "${meta.alias}"
    cpus params.threads
    // cannot use maxRetries based on exitcodes 137
    // due to the wf fail at the samtools step
    memory {12.GB * task.attempt}
    maxRetries 1
    errorStrategy 'retry'
    input:
        tuple val(meta), path(concat_seqs), path(fastcat_stats)
        path host_reference
        val common_minimap2_opts_list
    output:
        tuple(
            val(meta),
            path("*.unmapped.fastq.gz"),
            path("stats_unmapped"),
            env(n_seqs_passed_host_depletion),
            emit: fastq)
        tuple(
            val(meta),
            path("*.host.bam"),
            path("*.host.bam.bai"),
            emit: host_bam, optional:true)
        tuple(
            val(meta),
            path("*.unmapped.bam"),
            path("*.unmapped.bam.bai"),
            emit: no_host_bam, optional:true)
    script:
        def sample_id = "${meta.alias}"
        String fastcat_stats_outdir = "stats_unmapped"

        // use the file size to determine if minimap2 split-prefix is required. this is slightly conservative as newlines and headers will count against the reference size.
        String common_minimap2_opts = common_minimap2_opts_list.join(" ")
        String opt_split_prefix = host_reference.size() > 4e9  ? '--split-prefix tmp' : ""
    // Map reads against the host reference and take the unmapped reads for further analysis
    """
    minimap2 -t $task.cpus ${common_minimap2_opts} ${opt_split_prefix} -m 50 --secondary=no "${host_reference}" $concat_seqs \
    | samtools sort -@ ${task.cpus - 1}  --write-index -o "${sample_id}.all.bam##idx##${sample_id}.all.bam.bai" -
    # get unmapped reads & convert bam to fastq again
    samtools view -@ ${task.cpus - 1}  -b -f 4 --write-index \
        -o "${sample_id}.unmapped.bam##idx##${sample_id}.unmapped.bam.bai" \
        --unoutput "${sample_id}.host.bam##idx##${sample_id}.host.bam.bai" "${sample_id}.all.bam"
    samtools fastq -@ ${task.cpus - 1} -T '*' "${sample_id}.unmapped.bam" | bgzip -@ $task.cpus > "${sample_id}.unmapped.fastq.gz"

    # run fastcat on selected reads
    mkdir $fastcat_stats_outdir
    fastcat \
        -s "${sample_id}" \
        -f "${fastcat_stats_outdir}/per-file-stats.tsv" \
        --histograms histograms \
        "${sample_id}.unmapped.fastq.gz" > /dev/null
    # get number of reads after host removal
    n_seqs_passed_host_depletion=\$(awk 'NR==1{for (i=1; i<=NF; i++) {ix[\$i] = i}} NR>1 {c+=\$ix["n_seqs"]} END{print c}' \
        "${fastcat_stats_outdir}/per-file-stats.tsv")
    mv histograms/* $fastcat_stats_outdir
    """
}


workflow run_common {
    take:
        samples
        host_reference
        common_minimap2_opts
    main:
        reads = exclude_host_reads(samples,
            host_reference,
            common_minimap2_opts
        )
        samples_passed = reads.fastq.map {
            meta, fastq, stats, n_seqs_passed_host_depletion->
            [meta + [n_seqs_passed_host_depletion: n_seqs_passed_host_depletion as Integer], fastq, stats]
        }
        // Discard empty samples after host depletion
        branched = samples_passed
        | branch { meta, seqs, stats ->
            pass: meta.n_seqs_passed_host_depletion > 0
            fail: true
        }
        // log info about the failing samples
        branched.fail
            | subscribe {
                meta, seqs, stats ->
                def valid = meta['n_seqs_passed_host_depletion'] > 0
                if (!valid) {
                    log.warn "Found empty file after host depletion for sample: '${meta["alias"]}'."
                }
            onComplete: {
                log.warn "Empty files or those files whose reads have been discarded after host depletion " +
                "will not appear in the report and will be excluded from subsequent analysis."
            }
        }
        // Save the passing samples
        samples = branched.pass
    emit:
        samples
}