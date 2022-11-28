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
    storeDir "${params.store_dir}" 
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
    mv "${kmer_distribution}" database_dir
    """
}


process unpackTaxonomy {
    label "wfmetagenomics"
    cpus 1
    input:
        path taxonomy
    output:
        path "taxonomy_dir"
    """
    if [[ "${taxonomy}" == *.tar.gz ]]
    then
        mkdir taxonomy_dir
        tar xf "${taxonomy}" -C taxonomy_dir
    elif [ -d "${taxonomy}" ]
    then
        mv "${taxonomy}" taxonomy_dir
    else
        echo "Error: taxonomy is neither .tar.gz nor a dir"
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


process filterFastq {
    label "wfmetagenomics"
    maxForks params.threads  // no point having more inputs than processing threads
    cpus 2
    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("${sample_id}.${task.index}.fastq.gz"), emit: filtered
        tuple val(sample_id), path("${sample_id}.${task.index}.json"), emit: stats
    script:
        def max_length = "${params.max_len}"== null ? "-b ${params.max_len}" : ""
    shell:
    """
    fastcat \
        -a "${params.min_len}" \
        ${max_length} \
        -s "${sample_id}" \
        -r "${sample_id}.${task.index}.stats" \
        "${fastq}" | bgzip > "${sample_id}.${task.index}.fastq.gz"
    fastcat_histogram.py \
        --sample_id "${sample_id}" \
        "${sample_id}.${task.index}.stats" "${sample_id}.${task.index}.json"
    """
}


// scan step for accumulating fastcat stats
process progressiveStats {
    label "wfmetagenomics"
    maxForks 1
    cpus 1
    input: 
        path fastcat_stats
    output:
        path("stats.${task.index}")
    script: 
        if (fastcat_stats instanceof BlankSeparatedList) {
            new_input = fastcat_stats.getAt(0);
            state = fastcat_stats.getAt(-1)
        } else { 
            new_input = fastcat_stats;
            state = "stats.${task.index}"
        }
        output = "fastcat_stats.csv"
    """
    if [[ "${task.index}" == "1" ]]; then
        mkdir "$state";
    fi
    if [[ ! -f "$state/all_stats.json" ]]; then
        touch "$state/all_stats.json";
    fi
    add_jsons.py --new_file "${new_input}" --state "${state}/all_stats.json"
    mv "all_stats.json" "$state/all_stats.json"
    if [[ "${task.index}" != "1" ]]; then
         mv "$state" "stats.${task.index}"
    fi
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
    with open("$json/all_stats.json") as json_file:
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


process kraken2_client {
    label "wfmetagenomics"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    // retry if server responds out of resource
    errorStrategy = { task.exitStatus in [8] ? 'retry' : 'finish' }
    maxForks kraken_compute
    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path(reads), path("${sample_id}.kraken2.assignments.tsv"), path("${sample_id}.kraken2.report.txt"), emit: assignments
    script:
    """
    kraken2_client --port $params.port --report report.txt --sequence "${reads}" > "${sample_id}.kraken2.assignments.tsv"
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


// Combine kraken reports scan step
// Named a directory with a unique id to avoid dir name clash, but overwrite 'state' files in directory.
// Sample_id has to come from filename because can't scan tuples
// - maxForks is one to keep things in order, but this step can be slow.
process progressive_kreports {
    label "wfmetagenomics"
    maxForks 1 
    input:
        path kreport
        path taxonomy
    output:
        path("reports.${task.index}")
        val ("${sample_id}")
    script:
        if (kreport instanceof BlankSeparatedList){
            new_input =  kreport.getAt(0); state = kreport.getAt(-1)
        }
        else {
            new_input = kreport; state = "reports.${task.index}"
        }
        sample_id = "${new_input}".split(/\./)[0] 
    """
    if [[ "${task.index}" == "1" ]]; then
        mkdir "$state";
    fi
    if [[ ! -f "kreports/${sample_id}.kreport.txt" ]]; then
        touch "$state/${sample_id}.kreport.txt";
    fi
    combine_kreports_modified.py -r "$new_input"  "$state/${sample_id}.kreport.txt" -o ${sample_id}.kreport.txt --only-combined --no-headers
    mv ${sample_id}.kreport.txt $state/${sample_id}.kreport.txt
    if [[ "${task.index}" != "1" ]]; then
         mv "$state" "reports.${task.index}"
    fi
    """
}

// TODO: What is this process doing?!?!?
//       We end up remitting the kraken report directory and that gets used to build
//       our HTML report? Why isn't the .kreport_bracken_species.txt file 
//       always produced? And why do we not seem to care? How on Earth is it that the
//       .json files produced here end up in the progressive_kreports.out[0] directories
//       that get published? (Answer to the last one is they sometimes are, they sometimes
//       aren't because we've lost encapsulation by copying things through symlinks to
//       the input directory). 
process bracken {
    label "wfmetagenomics"
    cpus 1
    input:
        path kraken2_report
        val sample_id
        path database
        path taxonomy
    output:
        tuple val("${sample_id}"), path("*.kreport_bracken_species.txt"), emit: bracken_report, optional: true
        path kraken2_report, emit: report
    script:
    def awktab="awk -F '\t' -v OFS='\t'"
    """
    run_bracken.py \
        "${database}" \
        "${kraken2_report}/${sample_id}.kreport.txt" \
        "${params.bracken_length}" \
        "${params.bracken_level}" \
        "${sample_id}.bracken_report.txt"
    mv "${kraken2_report}/${sample_id}.kreport_bracken_species.txt" . || echo "no bracken report"

    ${awktab} '{ print \$2,\$6 }' "${sample_id}.bracken_report.txt" \
        | ${awktab} 'NR!=1 {print}' \
        | tee taxacounts.txt \
        | ${awktab} '{ print \$1 }' > taxa.txt
    taxonkit \
        --data-dir $taxonomy \
        lineage -R taxa.txt  > lineages.txt

    aggregate_lineages_bracken.py \
        -i "lineages.txt" -b "taxacounts.txt" \
        -u "${kraken2_report}/${sample_id}.kreport.txt" \
        -p "${sample_id}.kraken2"

    file1=`cat *.json`
    echo "{"'"$sample_id"'": "\$file1"}" >> "$sample_id.${task.index}.json"
    cp "${sample_id}.${task.index}.json" "${kraken2_report}/${sample_id}.json"
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
    script:
        report_name = "wf-metagenomics-report.html"
    """
    report.py \
        "${report_name}" \
        --versions versions \
        --params params.json \
        --summaries ${stats}/* \
        --lineages "${lineages}" \
        --vistempl template.html \
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
            batch_items = rebatchFastq(batch_items.transpose()).transpose()
        }

        // filter reads and calculate stats
        reads = filterFastq(batch_items)

        // progressive stats -- scan doesn't like tuple :/
        stats = progressiveStats.scan(
            reads.stats.map{ it -> it[1] })
   
        //  Stop file to input folder when read_limit stop condition is met.
        if (params.watch_path && params.run_indefinitely == false){
            stopCondition(stats)
        }

        // Run Kraken2
        kraken2_response = kraken2_client(reads.filtered)
    
        combined_kreport = progressive_kreports.scan(
            kraken2_response.map{ it -> it[3] },
            taxonomy)
        // TODO: what do we actually want back from this? See notes above `process bracken {}`
        br = bracken(combined_kreport, database, taxonomy)
        
        // report step
        //   Nasty: we assume br.report and stats are similarly ordered, they aren't.
        //   This doesn't matter too much as stats is read lengths and qualities: its
        //   not going to lead to grossly contradictory results. .merge() is deprecated
        //   so just ram the channels into the process. The alternative is needing to add
        //   keys into br.report and stats, but thats difficult because process.scan doesn't
        //   accept scan, because, ya know, reasons.
        versions = getVersions().collect()
        parameters = getParams().collect()
        stuff = stats
            .combine(versions)
            .combine(parameters)
            .combine(Channel.of(template))
        report = makeReport(br.report, stuff)
        
        // Stop server when all are processed
        stop_kraken_server(kraken2_response.collect())

        // output updating files as part of this pipeline
        // note, we don't overwrite the progressive kreport for now
        // TODO: see notes above `process bracken {}`
        output(report.report_html.mix(
            versions, parameters,
            progressive_kreports.out[0]),
        )

    emit:
        report.report_html  // just emit something
}
