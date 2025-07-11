import groovy.json.JsonBuilder


include { getParams } from '../../lib/common'

process abricateVersion {
    label "amr"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt", overwrite: true
    cpus 1
    memory "2 GB"
    input:
        path "input_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    cat input_versions.txt >> versions.txt
    abricate --version | sed 's/ /,/' >> "versions.txt"
    """
}

process getVersions {
    label "wfmetagenomics"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "versions.txt"
    cpus 1
    memory "2 GB"
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
    errorStrategy = 'retry'
    input:
        tuple val(meta), path(concat_seqs), path(fastcat_stats)
        path host_reference
        val common_minimap2_opts
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
        // use the file size to determine if minimap2 split-prefix is required. this is slightly conservative as newlines and headers will count against the reference size.
        def common_minimap2_opts = (host_reference.size() > 4e9 ) ? common_minimap2_opts + ['--split-prefix tmp'] : common_minimap2_opts
        common_minimap2_opts = common_minimap2_opts.join(" ")
        String fastcat_stats_outdir = "stats_unmapped"
        def per_read_stats = params.real_time ? "-r >(bgzip -c > $fastcat_stats_outdir/per-read-stats.tsv.gz)" : ""
    // Map reads against the host reference and take the unmapped reads for further analysis
    """
    minimap2 -t $task.cpus ${common_minimap2_opts} -m 50 --secondary=no "${host_reference}" $concat_seqs \
    | samtools sort -@ ${task.cpus - 1}  --write-index -o "${sample_id}.all.bam##idx##${sample_id}.all.bam.bai" -
    # get unmapped reads & convert bam to fastq again
    samtools view -@ ${task.cpus - 1}  -b -f 4 --write-index -o "${sample_id}.unmapped.bam##idx##${sample_id}.unmapped.bam.bai" --unoutput "${sample_id}.host.bam##idx##${sample_id}.host.bam.bai" "${sample_id}.all.bam"
        samtools fastq -@ ${task.cpus - 1} -T '*' "${sample_id}.unmapped.bam" | bgzip -@ $task.cpus > "${sample_id}.unmapped.fastq.gz"

    # run fastcat on selected reads
    mkdir $fastcat_stats_outdir
    fastcat \
        -s "${sample_id}" \
        -f $fastcat_stats_outdir/per-file-stats.tsv \
        --histograms histograms \
        ${per_read_stats} \
        "${sample_id}.unmapped.fastq.gz" > /dev/null
    # get number of reads after host removal
    n_seqs_passed_host_depletion=\$(awk 'NR==1{for (i=1; i<=NF; i++) {ix[\$i] = i}} NR>1 {c+=\$ix["n_seqs"]} END{print c}' \
        $fastcat_stats_outdir/per-file-stats.tsv)
    mv histograms/* $fastcat_stats_outdir
    """
}


// Process to collapse lineages info into abundance dataframes.
process createAbundanceTables {
    label "wfmetagenomics"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "abundance_table_*.tsv"
    cpus 1
    memory "2 GB"
    input:
        // lineages is a folder in the kraken2, but is a list of files in the minimap2 approach
        path "lineages/*"
        val taxonomic_rank
        val pipeline
    output:
        path("abundance_table_*.tsv"), emit: abundance_tsv

    """
    workflow-glue abundance_tables \
        --lineages lineages \
        --taxonomic_rank "${taxonomic_rank}" \
        --pipeline "${pipeline}"
    """
}

/* Extract reads in FASTQ from a list of IDs.
Use for example to output the unclassified reads.
 */
process publishReads {
    label "wfmetagenomics"
    publishDir "${params.out_dir}/${output_name}", mode: 'copy', pattern: "*.${output_name}.fq.gz", enabled: params.output_unclassified
    tag "${meta.alias}"
    cpus 1
    memory 4.GB
    input:
        tuple val(meta), path("reads.fq.gz"), path("ids.txt")
        val output_name
    output:
        path "${meta.alias}.${output_name}.fq.gz"
    script:
        """
        seqkit grep --pattern-file ids.txt reads.fq.gz -o "${meta.alias}.${output_name}.fq.gz"
        """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process publish {
    // publish inputs to output directory
    label "wfmetagenomics"
    cpus 1
    memory "2 GB"
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


workflow run_common {
    take:
        samples
        common_minimap2_opts
    main:
        common_versions = getVersions()
        if (params.amr){
            versions = abricateVersion(common_versions)
        } else{
            versions = common_versions
        }
        parameters = getParams()
        if (params.exclude_host){
            host_reference = file(params.exclude_host, checkIfExists: true)
            reads = exclude_host_reads(samples,
                host_reference,
                common_minimap2_opts
            )
            samples_host_depleted = reads.fastq.map {
                meta, fastq, stats, n_seqs_passed_host_depletion->
                [meta + [n_seqs_passed_host_depletion: n_seqs_passed_host_depletion as Integer], fastq, stats]
            }
            // Discard empty samples after host depletion
            branched = samples_host_depleted
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
            samples_passed = branched.pass
        } else{
            samples_passed = samples
        }
    emit:
        software_versions = versions
        parameters = parameters
        samples = samples_passed
}
