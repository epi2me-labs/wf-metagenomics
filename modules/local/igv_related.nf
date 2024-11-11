import groovy.json.JsonBuilder

// Process all samples together.
// Some references can be in just one sample but not in others
// We need the union because there is a common reference
process filter_references {
    label "wfmetagenomics"
    publishDir "${params.out_dir}/igv_reference", mode: 'copy', pattern: "reduced_reference.*"
    cpus params.threads
    memory "7 GB"
    input:
        path reference
        path "bam_flagstats/*"
    output:
        tuple(
            path("reduced_reference.fasta.gz"),
            path("reduced_reference.fasta.gz.fai"),
            path("reduced_reference.fasta.gz.gzi")
        )

    script:
    """
    # get found references for all the samples
    awk 'FNR > 1 && \$3 > ${params.abundance_threshold} && \$1 != "*" {print \$1}' bam_flagstats/*tsv \
    | sort -u \
    | seqkit grep -i -f - "${reference}" \
    | bgzip -@ $task.cpus > reduced_reference.fasta.gz
    samtools faidx reduced_reference.fasta.gz
    """
}