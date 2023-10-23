import groovy.json.JsonBuilder

process prepareSILVA {
    storeDir "${params.store_dir}/${params.database_set}"
    label "wfmetagenomics"
    cpus 2
    output:
        path "database/silva.fna", emit: reference
        path "taxonomy", emit: taxonomy
        path "database", emit: database
        path "seqid2taxid.map", emit: ref2taxid
        path "database1000mers.kmer_distrib", emit: bracken_dist
    script:
    """
    kraken2-build --db $params.database_set --special silva
    bracken-build -d $params.database_set -t "${task.cpus}" -l 1000
    BRACKEN_PATH=\$(dirname \$(which bracken))
    python \$BRACKEN_PATH/src/generate_kmer_distribution.py -i $params.database_set/database1000mers.kraken -o database1000mers.kmer_distrib
    # Move all the files following other default databases:
    mkdir database
    mv $params.database_set/*.k2d database/
    mv $params.database_set/library/silva.fna database/
    mv $params.database_set/seqid2taxid.map .
    mv $params.database_set/taxonomy taxonomy
    """
}

process abricateVersion {
    label "amr"
    cpus 1
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

process exclude_host_reads {
    label "wfmetagenomics"
    cpus params.threads
    input:
        tuple val(meta), path(concat_seqs), path(fastcat_stats)
        path host_reference
    output:
        tuple(
            val(meta),
            path("*.filtered.unmapped.fastq"),
            path("fastcat_stats_filtered"),
            emit: fastq)
        tuple(
            val(meta),
            path("*.host.bam"),
            path("*.host.bam.bai"),
            emit: host_bam)
    script:
        def sample_id = "${meta.alias}"
        def split = params.split_prefix ? '--split-prefix tmp' : ''
        String fastcat_stats_outdir = "fastcat_stats_filtered"
    // Map reads agains the host reference and take the unmapped reads for further analysis
    """
    minimap2 -t "${task.cpus}" ${split} -ax map-ont -m 50 --secondary no "${host_reference}" "${concat_seqs}" \
    | samtools view -h -b - | samtools sort -o "${sample_id}.all.bam" -
    samtools index "${sample_id}.all.bam"

    # get unmapped reads & convert bam to fastq again
    samtools view -b -f 4 "${sample_id}.all.bam" \
    | samtools fastq - > "${sample_id}.unmapped.fastq"

    # return host reads bam
    samtools view -b -F 4 "${sample_id}.all.bam" \
    | samtools sort -o "${sample_id}.host.bam" -
    samtools index "${sample_id}.host.bam"

    # run fastcat on selected reads
    mkdir $fastcat_stats_outdir
    fastcat \
            -s ${meta["alias"]} \
            -r $fastcat_stats_outdir/per-read-stats.tsv \
            -f $fastcat_stats_outdir/per-file-stats.tsv \
            "${sample_id}.unmapped.fastq" > "${sample_id}.filtered.unmapped.fastq"
    bgzip $fastcat_stats_outdir/per-read-stats.tsv 
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_host {
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
    echo "Writing host files"
    """
}

// Process to collapse lineages info into abundance dataframes.
process createAbundanceTables {
    label "wfmetagenomics"
    input:
        // lineages is a folder in the kraken2, but is a list of files in the minimap2 approach
        path "lineages/*"
        val taxonomic_rank
    output:
        path("abundance_table_*.tsv"), emit: abundance_tsv

    """
    workflow-glue abundance_tables \
        --lineages lineages \
        --taxonomic_rank "${taxonomic_rank}" \
        --pipeline "${params.classifier}"
    """
}


workflow run_common {
    take:
        samples
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
            reads = exclude_host_reads(samples, host_reference)
            samples = reads.fastq
            ch_to_publish = Channel.empty()
            ch_to_publish = ch_to_publish | mix (
            reads.host_bam | map { meta, bam, bai -> [bam, "host_bam"]},
            reads.host_bam | map { meta, bam, bai -> [bai, "host_bam"]},
            )
            ch_to_publish | output_host
        } else{
            samples
        }
    emit:
        software_versions = versions
        parameters = parameters
        samples
}
