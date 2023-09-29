import nextflow.util.BlankSeparatedList

include { run_amr } from '../modules/local/amr'
include { run_common } from '../modules/local/common'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process unpackDatabase {
    label "wfmetagenomics"
    cpus 1
    storeDir "${params.store_dir}/${database.simpleName}"
    input:
        path database
        path kmer_distribution
    output:
        path "database_dir"
    """
    if [[ "${database}" == *.tar.gz ]]
    then
        mkdir database_dir
        tar xf "${database}" -C database_dir
    elif [ -d "${database}" ]
    then
        mv "${database}" database_dir
    else
        echo "Error: database is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    cp "${kmer_distribution}" database_dir
    """
}


process determine_bracken_length {
    label "wfmetagenomics"
    input:
        path database
    output:
        env BRACKEN_LENGTH
    """
    if [[ -f "${database}"/database${params.bracken_length}mers.kmer_distrib ]]; then
        BRACKEN_LENGTH="${params.bracken_length}"
    else
        cd "${database}"
        BRACKEN_LENGTH=\$(ls -v1 *.kmer_distrib | tail -1 | sed -e "s@^database@@" -e "s@mers.kmer_distrib@@")
        cd ..
    fi
    """
}


process unpackTaxonomy {
    label "wfmetagenomics"
    cpus 1
    storeDir "${params.store_dir}/${taxonomy.simpleName}"
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
        unzip  "${taxonomy}" -d taxonomy_dir
    elif [ -d "${taxonomy}" ]
    then
        mv "${taxonomy}" taxonomy_dir
    else
        echo "Error: taxonomy is neither .tar.gz, .zip nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
} 


// Rebundle fastqs (this is mostly in case we're given one big file)
process rebatchFastq {
    label "wfmetagenomics"
    maxForks params.threads  // no point having more inputs than processing threads
    cpus 3
    input:
        tuple val(meta), path(fastq), path(stats)
    output:
        tuple(
            val(meta),
            path("fastq/${meta.alias}.part_*.fastq.gz"),
            path("fastcat_stats/${meta.alias}_part_*.tsv")
            )
    script:
        def sample_id = "${meta.alias}"
    """
    # Batch fasta file
    seqkit split2 "${fastq}" \
        -j ${task.cpus - 1} -s ${params.batch_size} \
        -e .gz -o "${sample_id}" \
        -O fastq
    # Batch stats file
    mkdir -p fastcat_stats
    tail -n +2 "${stats}/per-read-stats.tsv" | split -l "${params.batch_size}" \
    -d --additional-suffix=.tsv \
    --filter='sh -c "{ head -n1 "${stats}/per-read-stats.tsv"; cat; } > fastcat_stats/\$FILE"' - "${sample_id}"_part_
    """
}


// watch path stop condition, if params.read_limit is met will inject a stop file in to input folder.
process stopCondition { 
    label "wfmetagenomics"
    cpus 1 
    publishDir params.fastq, mode: 'copy', pattern: "*"
    input:
        path json
    output:
        path "STOP.fastq.gz", optional: true, emit: stop
    script:
        int threshold = params.read_limit
    """    
    #!/usr/bin/env python
    import json
    with open("$json") as json_file:
        state = json.load(json_file)
        total = 0 
        for k,v in state.items():
            total += v["total_reads"]
        if total >= $threshold:
            with open("STOP.fastq.gz", "x") as f:
                pass
    """
}


// Notes on CPU resource of kraken server and client:
//
// - server will use as much resource as number of clients running,
//   plus one extra thread => set maxForks of clients to kraken_clients - 1.
//   (not doing so gives gRPC server: "Server Threadpool Exhausted")
// - we need potentially one extra request thread to allow stop request
//   to be handled
// - cannot start so many clients (or other processes) such that
//   server never starts from Nextflow executor limit
// - we'd like to leave some resource for downstream processes such that we
//   get reporting updated frequently
// - this might all be considered a bit inefficient - but we are set
//   up for real-time dynamism not speed
//

kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process kraken_server {
    label "wfmetagenomics"
    cpus params.threads
    // short term solution: the following db require at least 8GB of memory,
    // normally in laptops with 16GB, docker only uses 8GB, so the wf stalls
    errorStrategy = {
        task.exitStatus == 137 & params.database_set in [
            'PlusPF-8', 'PlusPFP-8']? log.error("Error 137 while running kraken2_server, this may indicate the process ran out of memory.\nIf you are using Docker you should check the amount of RAM allocated to your Docker server.") : ''
        }
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        path database
    output:
        val true
    script:
        def memory_mapping = params.kraken2_memory_mapping ? '--memory-mapping' : ''
    """
    # we add one to requests to allow for stop signal
    kraken2_server \
        --port ${params.port} --host-ip ${params.host} \
        --max-requests ${kraken_compute + 1} --thread-pool ${params.server_threads}\
        --confidence ${params.kraken2_confidence}\
        --db ./${database} ${memory_mapping}
    """
}


// Filter reads, calculate some stats, and run kraken2 classification
//
// The per file statistics are calculated here rather than in a separate process
// in order to more easily define the parallelism (otherwise an independent stats
// job may queue up hundreds of jobs before kraken jobs, leaving the server consuming
// resource but not actually doing anything).
process kraken2_client {
    label "wfmetagenomics"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // retry if server responds out of resource
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    maxForks kraken_compute
    input:
        tuple val(meta), path(fastq), path(stats)
    output:
        tuple val("${meta.alias}"), path("${meta.alias}.kraken2.report.txt"), path("${meta.alias}.${task.index}.json"), path("kraken2.assignments.tsv")
    script:
        def sample_id = "${meta.alias}"
    """
    if [[ -f $stats ]]
    then
        stats_file="${stats}"
    else
        stats_file="${stats}/per-read-stats.tsv"
    fi

    workflow-glue fastcat_histogram \
        --sample_id "${sample_id}" \
        \$stats_file "${sample_id}.${task.index}.json"

    kraken2_client \
        --port ${params.port} --host-ip ${params.host} \
        --report report.txt \
        --sequence $fastq > "kraken2.assignments.tsv"
    tail -n +1 report.txt > "${sample_id}.kraken2.report.txt"
    """
}


process stop_kraken_server {
    label "wfmetagenomics"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // this shouldn't happen, but we'll keep retrying
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    input:
        val stop
    """
    kraken2_client --port $params.port --shutdown
    """
}


// Scan step for accumulating fastcat stats
//
// Nextflow scan does a silly thing where it feeds back the growing list of
// historical outputs. We only ever need the most recent output (the "state").
process progressive_stats {
    label "wfmetagenomics"
    maxForks 1
    cpus 1
    input: 
        path fastcat_stats
    output:
        path("all_stats.${task.index}")
    script:
        new_input = fastcat_stats instanceof BlankSeparatedList ? fastcat_stats.first() : fastcat_stats
        state = fastcat_stats instanceof BlankSeparatedList ? fastcat_stats.last() : "NOSTATE"
        output = "all_stats.${task.index}"
    """
    touch "${state}"
    workflow-glue add_jsons "${new_input}" "${state}" "${output}"
    """
}


// Combine kraken reports scan step
//
// Our state here is a directory of aggregated kraken2 reports, one file per
// sample. The state is named with the task.index so its unique, and also
// with the sample_id just analysed.
//
// Every new input to the process is a kraken report for a single sample. We
// don't want to update all aggregated reports, just the report for the sample
// of the latest report received. Scan cannot handle a tuple as a piece of
// state because the input channel is bodged to contain a list of the newest
// input and the growing list -- so we have to grotesquely pass the sample_id
// in a second channel (that is also accumulating).
process progressive_kraken_reports {
    label "wfmetagenomics"
    maxForks 1 
    publishDir path: "${params.out_dir}", mode: 'copy', pattern: "${new_state}", saveAs: {name -> "kraken"}, overwrite: true
    input:
        path kreport
        val sample_ids
        path taxonomy
    output:
        path("kraken.${task.index}.${sample_id}"), emit: reports
        val(sample_id), emit: sample_id
        path taxonomy
    script:
        def new_input = kreport instanceof List ? kreport.first() : kreport
        def state = kreport instanceof List ? kreport.last() : "NOSTATE"
        sample_id = sample_ids instanceof List ? sample_ids.first() : sample_ids
        new_state = "kraken.${task.index}.${sample_id}"
        // n.b where this is used below the files will have been moved, hence new_state
        old_input = "${new_state}/${sample_id}.kreport.txt"
    """
    if [[ "${task.index}" == "1" ]]; then
        mkdir "${state}"
    fi

    cp -r "${state}" "${new_state}" 
    touch "${old_input}"

    workflow-glue combine_kreports_modified \
        -r "${new_input}" "${old_input}" \
        -o "${sample_id}.kreport.txt" --only-combined --no-headers
    mv "${sample_id}.kreport.txt" "${new_state}/${sample_id}.kreport.txt"
    """
}


// Calculate up-to-date bracken information for latest analysed sample.
// 
// The kraken.scan process gives a directory of aggregated kraken2 reports
// together with the sample_id which was most recently updated. Here we rerun bracken
// on that sample and output the results to an updating directory of per sample
// bracken results.
process progressive_bracken {
    label "wfmetagenomics"
    cpus 2
    maxForks 1
    publishDir path: "${params.out_dir}", mode: 'copy', pattern: "${new_state}", saveAs: {name -> "bracken"}, overwrite: true
    input:
        path(inputs)
        val(sample_ids)
        tuple path(database), path(taxonomy), val(bracken_length), val(taxonomic_rank)
    output:
        path("${new_state}"), emit: reports
        val(sample_id), emit: sample_id
        tuple path(database), path(taxonomy), val(bracken_length), val(taxonomic_rank)
    script:
        new_state = "bracken.${task.index}"
        def kreports = inputs instanceof List ? inputs.first() : inputs
        def state = inputs instanceof List ? inputs.last() : "NOSTATE"
        sample_id = sample_ids instanceof List ? sample_ids.first(): sample_ids
        def awktab="awk -F '\t' -v OFS='\t'"
    """
    # run bracken on the latest kreports, is this writing some outputs
    # alongside the inputs? seems at least {}.kreport_bracken_species.txt
    # is written alongside the input
    workflow-glue run_bracken \
        "${database}" \
        "${kreports}/${sample_id}.kreport.txt" \
        "${bracken_length}" \
        "${taxonomic_rank}" \
        "${sample_id}.bracken_report.txt"

    # do some stuff...
    ${awktab} '{ print \$2,\$6 }' "${sample_id}.bracken_report.txt" \
        | ${awktab} 'NR!=1 {print}' \
        | tee taxacounts.txt \
        | ${awktab} '{ print \$1 }' > taxa.txt
    taxonkit lineage \
        -j ${task.cpus} \
        --data-dir $taxonomy \
        -R taxa.txt  > lineages.txt
    workflow-glue aggregate_lineages_bracken \
        -i "lineages.txt" -b "taxacounts.txt" \
        -u "${kreports}/${sample_id}.kreport.txt" \
        -p "${sample_id}.kraken2" \
        -r "${taxonomic_rank}" \

    file1=`cat *.json`
    echo "{"'"$sample_id"'": "\$file1"}" >> "bracken.json"

    # collate the latest bracken outputs into state
    if [[ "${task.index}" != "1" ]]; then
        cp -r "${state}" "${new_state}"
    else
        # make fresh directory
        mkdir "${new_state}"
    fi;

    # first output here is just for end user
    mv "${kreports}/${sample_id}.kreport_bracken_species.txt" "${new_state}" || echo "No bracken report"
    mv "bracken.json" "${new_state}/${sample_id}.json"
    """
}

// Concatenate kraken reports per read
process concatAssignments {
    label "wfmetagenomics"
    input:
        tuple (
        val(sample_id),
        path("classifications/kraken2.assignments.tsv")
        )
        path taxonomy
    output:
        tuple(
            val(sample_id),
            path("*_lineages.kraken2.assignments.tsv"), 
            emit: kraken2_reads_assignments
            )
    script:
    def sample_id = "${sample_id}"
    """
    cat classifications/* > all.sample.assignments.tsv
    # Run taxonkit to give users a more informative table
    taxonkit reformat  -I 3  --data-dir "${taxonomy}" -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" -F all.sample.assignments.tsv > "${sample_id}_lineages.kraken2.assignments.tsv"
    """

}

process makeReport {
    label "wfmetagenomics"
    maxForks 1
    cpus 1
    input:
        path lineages
        tuple(path(stats), path("versions/*"), path("params.json"), val(taxonomic_rank))
        path amr
    output:
        path "*.html", emit: report_html
    script:
        String workflow_name = workflow.manifest.name.replace("epi2me-labs/","")
        String report_name = "${workflow_name}-report.html"
        String amr_arg = amr.name != "OPTIONAL_FILE" ? "--amr ${amr}" : ""
    """
    workflow-glue report \
        "${report_name}" \
        --workflow_name ${workflow_name} \
        --versions versions \
        --params params.json \
        --read_stats ${stats} \
        --lineages "${lineages}" \
        --taxonomic_rank "${taxonomic_rank}" \
        --abundance_threshold "${params.abundance_threshold}"\
        --n_taxa_barplot "${params.n_taxa_barplot}"\
        $amr_arg
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


workflow kraken_pipeline {
    take:
        samples
        taxonomy
        database
        kmer_distribution
        taxonomic_rank
    main:
        opt_file = file("$projectDir/data/OPTIONAL_FILE")
        taxonomy = unpackTaxonomy(taxonomy)
        
        database = unpackDatabase(database, kmer_distribution)
        bracken_length = determine_bracken_length(database)
        

        // do we want to run a kraken server ourselves? 
        if (!params.external_kraken2) {
            kraken_server(database)
        }

        // maybe split up large files
        // sample_ids = samples.map { meta, reads, stats -> meta.alias }
        batch_items = samples
        if (params.batch_size != 0) {
            batch_items = rebatchFastq(batch_items.transpose())
                .transpose()
        }
        // filter host reads
        common = run_common(batch_items)
        software_versions = common.software_versions
        parameters = common.parameters
        samples_filtered = common.samples
        // Run Kraken2
        kraken2_client(samples_filtered)
    
        // progressive stats -- scan doesn't like tuple :/
        stats = progressive_stats.scan(
            kraken2_client.out.map { id, report, json, assignments -> json },
        )

        kraken2_client.out.multiMap {
            id, rep, json, assignments ->
                sample_id: id
                report: rep
            }
            .set { scan_input }

        progressive_kraken_reports.scan(
            scan_input.report, scan_input.sample_id, 
            taxonomy)

        database
            .combine(taxonomy)
            .combine(bracken_length)
            .combine(Channel.of(taxonomic_rank))
            .first() // To ensure value channel for scan
            .set {bracken_inputs}

        progressive_bracken.scan(
            progressive_kraken_reports.out.reports, progressive_kraken_reports.out.sample_id,
            bracken_inputs)

        // Process AMR
        if (params.amr) {
            run_amr = run_amr(
                batch_items,
                "${params.amr_db}",
                "${params.amr_minid}",
                "${params.amr_mincov}"
            )
            amr_reports = run_amr.reports
        } else {
            // use first() to coerce this to a value channel
	        amr_reports = Channel.fromPath("$projectDir/data/OPTIONAL_FILE", checkIfExists: true).first()
	}

        // report step
        // Nasty: we assume br.report and stats are similarly
        // ordered. This is find because everything has been through scan and
        // is therefore necessarily ordered. This could be tidied up by passing
        // through a job id but that seems unneccessary at this point.
        stuff = stats
            .combine(software_versions)
            .combine(parameters)
            .combine(Channel.of(taxonomic_rank))
        
    
        report = makeReport(
            progressive_bracken.out.reports,
            stuff,
            amr_reports
        )
        
        // output updating files as part of this pipeline
        ch_to_publish = Channel.empty()
        | mix(
            software_versions,
            parameters,
            report.report_html,
        )
        | map { [it, null] }

        if (params.include_kraken2_assignments) {
            kraken2_assignments = concatAssignments(
                kraken2_client.out.map{ 
                    id, report, json, assignments -> tuple(id, assignments) 
                }.groupTuple(), taxonomy)
            ch_to_publish = ch_to_publish | mix (
            kraken2_assignments.kraken2_reads_assignments | map {
                 id, kraken_reads_classification -> [kraken_reads_classification, "kraken_reads_assignments"]},
            )
        }

        ch_to_publish | output


        // Stop server when all are processed
        if (!params.external_kraken2) {
            stop_kraken_server(kraken2_client.out.collect())
        }

        //  Stop file to input folder when read_limit stop condition is met.
        if (params.watch_path && params.read_limit){
            stopCondition(stats).first().subscribe {
            def stop = file(params.fastq).resolve("STOP.fastq.gz")
                log.info "Creating STOP file: '$stop'"
            }
        }

    emit:
        report.report_html  // just emit something
}

workflow.onComplete {
    def stop = file(params.fastq).resolve("STOP.fastq.gz")
    if (stop.exists()) {
        stop.delete()
        log.info "Deleted STOP file: '$stop'"
    }
}
