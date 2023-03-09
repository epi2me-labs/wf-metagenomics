import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList


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
        tuple val(sample_id), path(fastq)
    output:
        tuple(
            val(sample_id),
            path("fastq/${sample_id}.part_*.fastq.gz"))
    shell:
    """
    seqkit split2 ${fastq} \
        -j ${task.cpus - 1} -s ${params.batch_size} \
        -e .gz -o "${sample_id}" \
        -O fastq
    """
}


// watch path stop condition, if params.read_limit is met will inject a stop file in to input folder.
process stopCondition { 
    label "wfmetagenomics"
    cpus 1 
    publishDir "${params.fastq}/STOP/", mode: 'copy', pattern: "*"
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
//   plus one extra thread => set maxForks of clients to threads - 1.
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

kraken_compute = params.threads == 1 ? 1 : params.threads - 1

process kraken_server {
    errorStrategy 'ignore'
    label "wfmetagenomics"
    cpus params.threads
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        path database
    output:
        val true
    script:
    """
    # we add one to requests to allow for stop signal
    kraken2_server \
        --max-requests ${kraken_compute + 1} --port ${params.port} \
        --db ./${database}/
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
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("${sample_id}.kraken2.report.txt"), path("${sample_id}.${task.index}.json")
    script:
        def max_length = "${params.max_len}"== null ? "-b ${params.max_len}" : ""
    """

    fastcat \
        -a "${params.min_len}" \
        ${max_length} \
        -s "${sample_id}" \
        -r "${sample_id}.${task.index}.stats" \
        "${fastq}" > filtered.fastq
    workflow-glue fastcat_histogram \
        --sample_id "${sample_id}" \
        "${sample_id}.${task.index}.stats" "${sample_id}.${task.index}.json"

    kraken2_client \
        --port $params.port --report report.txt \
        --sequence filtered.fastq > "${sample_id}.kraken2.assignments.tsv"
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
        tuple path(database), path(taxonomy), val(bracken_length)
    output:
        path("${new_state}"), emit: reports
        val(sample_id), emit: sample_id
        tuple path(database), path(taxonomy), val(bracken_length)
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
        "${params.bracken_level}" \
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
        -p "${sample_id}.kraken2"

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


process makeReport {
    label "wfmetagenomics"
    maxForks 1
    cpus 1
    input:
        path lineages
        tuple(path(stats), path("versions/*"), path("params.json"), path("template.html"))
    output:
        path "wf-metagenomics-*.html", emit: report_html
        path "wf-metagenomics-counts.tsv", emit: counts_tsv
        path "wf-metagenomics-rarefied.tsv", emit: rarefied_tsv
    script:
        report_name = "wf-metagenomics-report.html"
    """
    workflow-glue report \
        "${report_name}" \
        --versions versions \
        --params params.json \
        --summaries ${stats} \
        --lineages "${lineages}" \
        --vistempl template.html \
        --rank "${params.bracken_level}"
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


workflow kraken_pipeline {
    take:
        taxonomy
        database
        kmer_distribution
        template
    main:
 
        taxonomy = unpackTaxonomy(taxonomy)
        
        database = unpackDatabase(database, kmer_distribution)
        bracken_length = determine_bracken_length(database)
        kraken_server(database)
        input = file("${params.fastq}")
        if (input.isFile()){
            if (params.watch_path) {
                throw new Exception("Watch path can only be used with input directories")
            }
            initial_items = channel.fromPath("${params.fastq}")
            batch_items = initial_items
                .map{ it -> tuple(it.simpleName, it) }
        }
        else if (input.isDirectory()) {
            log.info("")
            log.info("Input directory assumed to be containing one or more directories containing fastq files.")
            initial_items = Channel.fromPath("${params.fastq}/**/*f*q*")
            all_items = initial_items
            if (params.watch_path){
                added_items = Channel
                    .watchPath("${params.fastq}/**/*f*q*")
                    .until{ it -> it.name == 'STOP.fastq.gz' }
                all_items = initial_items.concat(added_items)
            }
            named = all_items
                .map{ it -> tuple("$it".split(/\//)[-2], it) }
            batch_items = named
            if (params.sample_sheet){
                // Use sample sheet to name samples and batch
                sample_sheet = get_sample_sheet(params.sample_sheet)
                batch_items = named
                    .combine(sample_sheet, by: [0])
                    .map { it -> tuple(it[2], it[1]) }
            }
        }
        else {
            throw new Exception("--fastq should be a file or directory.")
        }

        // maybe split up large files
        if (params.batch_size != 0) {
            batch_items = rebatchFastq(batch_items.transpose())
                .transpose()
        }

        // Run Kraken2
        kraken2_client(batch_items)
    
        // progressive stats -- scan doesn't like tuple :/
        stats = progressive_stats.scan(
            kraken2_client.out.map { id, report, json -> json },
        )

        kraken2_client.out.multiMap {
            id, rep, json ->
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
            .first() // To ensure value channel for scan
            .set {bracken_inputs}

        progressive_bracken.scan(
            progressive_kraken_reports.out.reports, progressive_kraken_reports.out.sample_id,
            bracken_inputs)

        // report step
        // Nasty: we assume br.report and stats are similarly
        // ordered. This is find because everything has been through scan and
        // is therefore necessarily ordered. This could be tidied up by passing
        // through a job id but that seems unneccessary at this point.
        versions = getVersions().collect()
        parameters = getParams().collect()
        stuff = stats
            .combine(versions)
            .combine(parameters)
            .combine(Channel.of(template))
        
        report = makeReport(progressive_bracken.out.reports, stuff)
        
        // output updating files as part of this pipeline
        output(report.report_html.mix(
            versions, parameters, report.counts_tsv, report.rarefied_tsv),
        )

        // Stop server when all are processed
        stop_kraken_server(kraken2_client.out.collect())

        //  Stop file to input folder when read_limit stop condition is met.
        if (params.watch_path && params.read_limit){
            stopCondition(stats)
        }

    emit:
        report.report_html  // just emit something
        report.rarefied_tsv
        report.counts_tsv
}
