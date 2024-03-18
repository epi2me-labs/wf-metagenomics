import nextflow.util.BlankSeparatedList

include { run_amr } from '../modules/local/amr'
include {
    run_common;
    createAbundanceTables;
    output as output_results;
} from '../modules/local/common'
OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


// Rebundle fastqs (this is mostly in case we're given one big file)
process rebatchFastq {
    label "wfmetagenomics"
    tag "${meta.alias}"
    maxForks params.threads  // no point having more inputs than processing threads
    cpus 2
    memory "2GB"
    input:
        tuple val(meta), path(concat_seqs), path(stats)
    output:
        tuple(
            val(meta),
            path("fastq/*.part_*.fastq.gz"),
            path("fastq/fastcat_stats/*.tsv.gz")
        )
    script:
        def sample_id = "${meta.alias}"
    """
    seqkit split2 $concat_seqs \
        -j ${task.cpus - 1} -s ${params.batch_size} \
        -e .gz -o "${sample_id}" \
        -O fastq
    # Batch stats file
    # run fastcat on each of the batch files:
    # don't need extra_args because the sequences have already passed fastcat during ingress
    mkdir -p fastq/fastcat_stats
    for f in \$(ls fastq);
        do
        if [[ fastq/\$f = *.fastq.gz ]];
        then
            batch_part=(\${f/.fastq.gz/})
            fastcat \
                -s \$batch_part \
                -r >(bgzip -c > "fastq/fastcat_stats/\$batch_part-per-read-stats.tsv.gz") \
                -f "fastq/fastcat_stats/\$batch_part-per-file-stats.tsv" \
                fastq/\$f --histograms \$batch_part-histograms
        fi
    done
    """
}

// watch path stop condition, if params.read_limit is met will inject a stop file in to input folder.
process stopCondition { 
    label "wfmetagenomics"
    cpus 1 
    memory "2GB"
    publishDir = [
        path: { params.fastq ? params.fastq: params.bam },
        mode: 'copy', pattern: "*"
    ]
    input:
        path json
    output:
        path "STOP.*", optional: true, emit: stop
    script:
        int threshold = params.read_limit
        String input_path = params.fastq ? "fastq": "" // just pass something so that python can evaluate True or False
    """    
    #!/usr/bin/env python
    import json
    with open("$json") as json_file:
        state = json.load(json_file)
        total = 0 
        for k,v in state.items():
            total += v["total_reads"]
        if total >= $threshold:
            if "$input_path":
                with open("STOP.fastq.gz", "x") as f:
                    pass
            else:
                with open("STOP.bam", "x") as f:
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
// - must make sure to not start so many clients (or other processes)
//   such that the server never starts due to the Nextflow executor limit
// - we'd like to leave some resource for downstream processes such that we
//   get reporting updated frequently
// - this might all be considered a bit inefficient - but we are set
//   up for real-time dynamism not speed
//

kraken_compute = params.kraken_clients == 1 ? 1 : params.kraken_clients - 1

process kraken_server {
    label "wfmetagenomics"
    memory {
        "${hash_size + 1e9} B " // leave buffer of 2GB, it's just the db
    }
    cpus {
        // if the task executor is local, ensure that we cannot start the
        // server with the same number of threads as the local executor limit.
        // as it will not be possible to run a client at the same time and
        // indefinitely block the workflow. in this case, we'll ensure that
        // we can always run at least one kraken client (the default).
        if (task.executor == "local") {
            // attempt to count local max executor and fall back to local cpus
            Integer max_local_threads = workflow.session.config?.executor?.$local?.cpus ?: \
                Runtime.getRuntime().availableProcessors()

            if (max_local_threads == 1) {
                throw new Exception("Cannot run kraken_server and kraken_client at the same time as the local executor appears to be configured with only one CPU.")
            }
            else if (max_local_threads == 2) {
                // run the server single threaded and expect one client
                log.info("Automatically set kraken2 classification server threads to 1 to ensure a classification client can be run.")
                1
            }
            else {
                // remove one thread for at least one client, and another for other business
                Integer server_threads = Math.min(params.server_threads, max_local_threads - 2)
                if (server_threads != params.server_threads) {
                    String msg = "Adjusted kraken2 classification server threads ${params.server_threads} -> ${server_threads} to ensure a classification client can be run."
                    if (params.kraken_clients > 2) {
                        // escalate to warn as fewer clients than desired are running
                        log.warn(msg + "\nYou requested ${params.kraken_clients} kraken clients, but only one can be run. Decrease the number of kraken2 server threads to allow more clients, decrease the number of desired kraken clients, or run this workflow on a device with more CPUs.")
                    }
                    else {
                        log.info(msg)
                    }
                }
                server_threads
            }
        }
        else {
            // if we're not local, just give the user what they want
            params.server_threads
        }
    }
    // short term solution: the following db require at least 8GB of memory,
    // normally in laptops with 16GB, docker only uses 8GB, so the wf stalls
    errorStrategy {
        task.exitStatus == 137 & params.database_set in [
            'Standard-8', 'PlusPF-8', 'PlusPFP-8']? log.error("Error 137 while running kraken2_server, this may indicate the process ran out of memory.\nIf you are using Docker you should check the amount of RAM allocated to your Docker server.") : ''
        }
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        path database
        val hash_size
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
    tag "${meta.alias}"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // retry if server responds out of resource
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    maxForks kraken_compute
    cpus 1
    memory "2GB"
    input:
        tuple val(meta), path(concat_seqs), path(stats)
    output:
        tuple val("${meta.alias}"), path("${meta.alias}.kraken2.report.txt"), path("${meta.alias}.${task.index}.json"), path("kraken2.assignments.tsv")
    script:
        def sample_id = "${meta.alias}"
    """
    if [[ -f $stats ]]
    then
        stats_file="${stats}" # batch case
    else
    # -L is needed because fastcat_stats is a symlink
        stats_file=\$(find -L "${stats}" -name "*.tsv.gz" -exec ls {} +)
    fi
    workflow-glue fastcat_histogram \
        --sample_id "${sample_id}" \
        \$stats_file "${sample_id}.${task.index}.json"

    kraken2_client \
        --port ${params.port} --host-ip ${params.host} \
        --report report.txt \
        --sequence $concat_seqs > "kraken2.assignments.tsv"
    tail -n +1 report.txt > "${sample_id}.kraken2.report.txt"
    """
}


process stop_kraken_server {
    label "wfmetagenomics"
    cpus 1
    memory "2GB"
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
    memory "2GB"
    input: 
        path stats
    output:
        path("all_stats.${task.index}")
    script:
        def new_input = stats instanceof BlankSeparatedList ? stats.first() : stats
        def state = stats instanceof BlankSeparatedList ? stats.last() : "NOSTATE"
        String output = "all_stats.${task.index}"
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
    tag "${sample_id}"
    maxForks 1 
    cpus 1
    memory "2GB"
    publishDir path: "${params.out_dir}", mode: 'copy', pattern: "${new_state}", saveAs: {name -> "kraken"}, overwrite: true
    input:
        path kreport
        val sample_ids
    output:
        path("kraken.${task.index}.${sample_id}"), emit: reports
        val(sample_id), emit: sample_id
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
    tag "${sample_id}"
    cpus 1
    memory "2GB"
    maxForks 1
    publishDir path: "${params.out_dir}", mode: 'copy', pattern: "${new_state}", saveAs: {name -> "bracken"}, overwrite: true
    input:
        path(inputs)
        val(sample_ids)
        tuple path(database), path(taxonomy), path(bracken_length_file), val(taxonomic_rank)
    output:
        path("${new_state}"), emit: reports
        val(sample_id), emit: sample_id
        tuple path(database), path(taxonomy), path(bracken_length_file), val(taxonomic_rank)
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
    BRACKEN_LENGTH=\$(cat "${bracken_length_file}")
    workflow-glue run_bracken \
        "${database}" \
        "${kreports}/${sample_id}.kreport.txt" \
        \$BRACKEN_LENGTH \
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

    file1=\$(find -name '*.json' -exec cat {} +)
    echo "{"'"$sample_id"'": \$file1}" >> "bracken.json"

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
    tag "${sample_id}"
    maxForks 1
    cpus 1
    memory "2GB"
    input:
        tuple (
            val(sample_id),
            path("classifications/?_kraken2.assignments.tsv")
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
    find -L classifications -name '*_kraken2.assignments.tsv' -exec cat {} + > all.sample.assignments.tsv
    # Run taxonkit to give users a more informative table
    taxonkit reformat  -I 3  --data-dir "${taxonomy}" -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" -F all.sample.assignments.tsv > "${sample_id}_lineages.kraken2.assignments.tsv"
    """

}

process makeReport {
    label "wf_common"
    maxForks 1
    cpus 1
    memory "2GB" //depends on the number of different species identified that tables may be bigger.
    input:
        path lineages
        path abundance_table
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
        --abundance_table "${abundance_table}" \
        --taxonomic_rank "${taxonomic_rank}" \
        --abundance_threshold "${params.abundance_threshold}"\
        --pipeline "real_time" \
        --n_taxa_barplot "${params.n_taxa_barplot}"\
        $amr_arg
    """
}


workflow real_time_pipeline {
    take:
        samples
        taxonomy
        database
        bracken_length
        taxonomic_rank
    main:
        opt_file = file("$projectDir/data/OPTIONAL_FILE")

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
        // do we want to run a kraken server ourselves?
        def database_main_file_size = database.resolve('hash.k2d').size()
        if (!params.external_kraken2) {
            kraken_server(database, database_main_file_size)
        }
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
            scan_input.report,
            scan_input.sample_id
        )
        
        database
            .combine(taxonomy)
            .combine(bracken_length)
            .combine(taxonomic_rank)
            .first() // To ensure value channel for scan
            .set {bracken_inputs}

        progressive_bracken.scan(
            progressive_kraken_reports.out.reports,
            progressive_kraken_reports.out.sample_id,
            bracken_inputs
        )

        // Process AMR
        if (params.amr) {
            run_amr = run_amr(
                samples_filtered,
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
        // ordered. This is fine because everything has been through scan and
        // is therefore necessarily ordered. This could be tidied up by passing
        // through a job id but that seems unnecessary at this point.
        basic_report_components = stats
            .combine(software_versions)
            .combine(parameters)
            .combine(taxonomic_rank)
        
        abundance_tables = createAbundanceTables(
            progressive_bracken.out.reports,
            taxonomic_rank, 'real_time')
        
        report = makeReport(
            progressive_bracken.out.reports,
            abundance_tables.abundance_tsv,
            basic_report_components,
            amr_reports
        )
        
        // output updating files as part of this pipeline
        ch_to_publish = Channel.empty()
        | mix(
            software_versions,
            parameters,
            report.report_html,
            abundance_tables.abundance_tsv
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

        ch_to_publish | output_results


        // Stop server when all are processed
        if (!params.external_kraken2) {
            stop_kraken_server(kraken2_client.out.collect())
        }

        //  Stop file to input folder when read_limit stop condition is met.
        if (params.real_time && params.read_limit){
            stopCondition(stats).first().subscribe {
                def stop = params.fastq ? file(params.fastq).resolve("STOP.fastq.gz") : \
                    file(params.bam).resolve("STOP.bam")
                log.info "Creating STOP file: '$stop'"
            }
        }

    emit:
        report.report_html  // just emit something
}

workflow.onComplete {
    // Every closure defined with workflow.onComplete will run at the end of the main
    // workflow (this is also true when it was defined below a sub-workflow in a
    // different `.nf` file). We thus need to check the classifier again, as `params.fastq`/`params.bam`
    // might be `null` in the minimap case.
    if ((params.classifier == "kraken2") && params.real_time) {
        def stop = params.fastq ? file(params.fastq).resolve("STOP.fastq.gz") : \
            file(params.bam).resolve("STOP.bam")
        if (stop.exists()) {
            stop.delete()
            log.info "Deleted STOP file: '$stop'"
        }
    }
}
