import groovy.json.JsonBuilder

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
    output:
        path("abundance_table_*.tsv"), emit: abundance_tsv

    """
    workflow-glue abundance_tables \
        --lineages lineages \
        --taxonomic_rank "${taxonomic_rank}" \
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


process makeReport {
    label "wf_common"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "${report_name}"
    cpus 1
    memory 4.GB
    input:
        val wf_version
        val metadata
        path(stats, stageAs:"stats/stats_*")
        path "abundance_table.tsv"
        path "alignment_stats/*"
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
        String align_stats = params.minimap2_by_reference ? "--align_stats alignment_stats" : ""
        String amr = params.amr as Boolean ? "--amr amr" : ""
    """
    echo '${metadata}' > metadata.json
    workflow-glue report \
        "${report_name}" \
        --workflow_name ${workflow_name} \
        --versions versions \
        --params params.json \
        --workflow_version $wf_version \
        --metadata metadata.json \
        --read_stats $stats \
        --lineages lineages \
        --abundance_table "abundance_table.tsv" \
        --taxonomic_rank "${taxonomic_rank}" \
        --abundance_threshold "${params.abundance_threshold}" \
        --n_taxa_barplot "${params.n_taxa_barplot}" \
        ${align_stats} \
        ${amr}
    """
}
