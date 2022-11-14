import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList


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


process combineFilterFastq {
    label "wfmetagenomics"
    forks = params.threads.intdiv(6)*5 
    maxForks forks <= 0 ? 1 : forks
    cpus 1
    input:
        tuple val(sample_id), path(directory)
        path database
    output:
        tuple(
            val(sample_id),
            path("*.fastq.gz"),
            emit: filtered)
        tuple(
            val(sample_id),
            path("*.json"),
            emit: stats)
    script:
        def listFiles = directory.join(" ")
    shell:
    """
    fastcat \
        -a "${params.min_len}" \
        -b "${params.max_len}" \
        -q 10 \
        -s "${sample_id}" \
        -r "${sample_id}.${task.index}.stats" \
        ${listFiles} | bgzip > "${sample_id}.${task.index}.fastq.gz"
    fastcat_histogram.py --file "${sample_id}.${task.index}.stats" --sample_id "$sample_id" 
    """
}


process extractKraken2Reads {
    label "wfmetagenomics"
    forks = params.threads.intdiv(6)*5 
    maxForks forks <= 0 ? 1 : forks
    cpus 1
    input:
        tuple val(sample_id), path(kraken_assignments), path(kraken_report), path(reads)
    output:
        tuple(
            val(sample_id),
            path("*extracted.fastq"),
            emit: extracted)
    script:
        def taxids = (params.kraken2filter as String).replaceAll(',',' ').replaceAll("'","")
        def policy = params.kraken2exclude ? '--exclude' : ''  
    """
    
    extract_kraken_reads.py \
        -k "${kraken_assignments}" \
        -r "${kraken_report}" \
        -s1 "${reads}" \
        -o "${sample_id}.kraken2.extracted.fastq" \
        -t ${taxids} \
        --fastq-output \
        --include-children
        ${policy}
    """
}


process bracken {
    label "wfmetagenomics"
    forks = params.threads.intdiv(6)*5 
    maxForks forks <= 0 ? 1 : forks
    cpus 1
    input:
        path kraken2_report
        val sample_id
        path database
        path taxonomy
    output:
        tuple val("${sample_id}"), path("*.kreport_bracken_species.txt"),emit: bracken_report, optional: true
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


process makeReport {
    label "wfmetagenomics"
    maxForks 1
    cpus 1
    input:
        path lineages
        path stats
        path "versions/*"
        path "params.json"
        path template
    output:
        path "wf-metagenomics-*.html", emit: report_html
    script:
        report_name = "wf-metagenomics-" + params.report_name + '.html'
    """
    report.py \
        "${report_name}" \
        --versions versions \
        --params params.json \
        --summaries ${stats}/* \
        --lineages "${lineages}" \
        --vistempl "${template}" \
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
    sample_id = "${new_input}".split(/\./)[0] 
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




// the server processes
process kraken_server {
    errorStrategy 'ignore'
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    label "wfmetagenomics"
    input:
        path database
    output:
        val true
    script:
    """
    cd $database
    kraken2_server --db ./ --port $params.port || exit 0
    """
}


process kraken2_client {
    errorStrategy 'retry'
    maxErrors 5
    label "wfmetagenomics"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    maxForks 1 
    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path(reads), path("*.tsv"), path("*kraken2_report.txt"), emit: assignments
    script:
    """
    kraken2_client --port $params.port --sequence "${reads}" > "${sample_id}.kraken2.assignments.tsv"
    kraken2_client --port $params.port --report "tmp.txt" --sequence "${reads}"
    tail -n +1 "tmp.txt" > "${sample_id}.kraken2_report.txt"
    """

}

// combine kraken reports scan step
// Named a directory with a unique id to avoid dir name clash, but overwrite 'state' files in directory.
// Sample_id has to come from filename because can't scan tuples
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
            else { new_input = kreport; state = "reports.${task.index}" }
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


process taxon_kit {
    label "wfmetagenomics"
    input:
        tuple val(sample_id), path(reads), path(tsv), path(report)
        path taxonomy
    output:
        tuple(
            val(sample_id),
            path(tsv),
            path("*lineages.txt"),
            path("*lineages.json"),
            path(report))
    script:
    // this step takes a lot of time because of the extract_kraken python script so would be good to just output these files from the kraken step
    // can also just remove these steps but no classified fastq file output. 
    """
    awk -F '\\t' '{print \$3}' "${tsv}" > taxids.tmp
    taxonkit \
        --data-dir "${taxonomy}" \
        lineage -R taxids.tmp \
        | aggregate_lineages.py -p "${sample_id}.kraken2"
    """
}


process stop_kraken_server {
    label "wfmetagenomics"
    containerOptions {workflow.profile != "singularity" ? "--network host" : ""}
    input:
        val stop
    """
    kraken2_client --port $params.port --shutdown
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


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output_dir {
    // publish inputs to output directory
    label "wfmetagenomics"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path ("*")
    """
    cp $fname/* . 
    echo "Writing output files"
    """
}


process mergeKrakenFiltered {
    label "wfmetagenomics"
    cpus 1
    input:
        tuple val(sample_id), path("files*.gz")
    output:
        tuple (val(sample_id), path("*.fastq"), emit: reads)
    """
    fastcat files* > "${sample_id}.kraken2.extracted.fastq"
    """
}


process catAssignmentsprogressive {
    label "wfmetagenomics"
    cpus 1
    input:
        path assignments
        path taxonomy
    output:
        path "assignments.${task.index}", emit: assignments
        path taxonomy
    script:
    if (assignments instanceof BlankSeparatedList){
            new_input =  assignments.getAt(0); state = assignments.getAt(-1)
            }
            else { new_input = assignments; state = "assignments.${task.index}" }
    sample_id = "${new_input}".split(/\./)[0]

    """
  
    if [[ "${task.index}" == "1" ]]; then
        mkdir "$state";
    fi
    if [[ ! -f "$state/${sample_id}.kraken2.tsv" ]]; then
        touch "$state/${sample_id}.kraken2.tsv";
    fi
    cat "$new_input" "$state/${sample_id}.kraken2.tsv" > "${sample_id}.kraken2.tsv"
    mv "${sample_id}.kraken2.tsv" "$state/${sample_id}.kraken2.tsv"
    awk -F '\\t' '{print \$3}' "$state/${sample_id}.kraken2.tsv" > "taxids.tmp"
    taxonkit \
        --data-dir $taxonomy \
        lineage -R taxids.tmp \
        | aggregate_lineages.py -p ${sample_id}.kraken2
    mv "${sample_id}.kraken2.lineages.txt" "$state/${sample_id}.kraken2.lineages.txt"
    if [[ "${task.index}" != "1" ]]; then
        mv $state assignments.${task.index}
    fi
    

    """
}


workflow kraken_pipeline {
    take:
        reference
        refindex
        ref2taxid
        taxonomy
        database
        kmer_distribution
        template
    main:
        
        outputs = []
        taxonomy = unpackTaxonomy(taxonomy)
        
        database = unpackDatabase(
                database,
                kmer_distribution
        )
        kraken_server(database)
        
        input = file("${params.fastq}")
        if (input.isFile() && params.watch_path){
            throw new Exception("Watch path can only be used with input directories")
        }
        if (input.isFile()) {
            initial_items = channel.fromPath("${params.fastq}")
            batch_items = initial_items.map{ it -> return tuple(it.simpleName, it) }
        }
        if (input.isDirectory()) {
            log.info("")
            log.info("Input directory assumed to be containing one or more directories containing fastq files.")
            initial_items = Channel.fromPath("${params.fastq}/**/*f*q*")
        }
        all_items = initial_items
        if (params.watch_path){
            added_items = Channel.watchPath("${params.fastq}/**/*f*q*").until{ file->file.name == 'STOP.fastq.gz' }
            all_items = initial_items.concat(added_items)
         
        }   

        if (input.isDirectory()) {
            if (params.sample_sheet){
            //check sample sheet
            sample_sheet = get_sample_sheet(params.sample_sheet)
            //Use sample sheet to name samples and batch
            batch_items = all_items.buffer( size:params.batch_size )
                                        .map{ it -> return tuple("$it".split(/\//)[-2], it) }
                                        .combine(sample_sheet, by: [0])
                                        .map { it -> tuple(it[2], it[1]) }
            }else{
                //Use name of directories as sample_names and batch
                batch_items = all_items.buffer( size:params.batch_size ).map{ it -> return tuple("$it".split(/\//)[-2], it) }
            }
        }
        // combine reads
        reads = combineFilterFastq(batch_items, database)

        // progressive stats 
        just_stats = reads.stats.map{ it -> return it[1]}
        stats = progressiveStats.scan(just_stats)
   
        //  Stop file to input folder when read_limit stop condition is met. Only used when --watch_path is true
        if (params.watch_path && params.run_indefinitely == false){
            stopCondition(stats)
        }

        // Run Kraken2
        kraken2_response = kraken2_client(reads.filtered)
    
        combined_kreport = progressive_kreports.scan(kraken2_response.map { it -> it[3]}, taxonomy)
        kr2 = taxon_kit(kraken2_response.assignments, taxonomy)
        br = bracken(
                combined_kreport,
                database, 
                taxonomy
            )
        outputs += [br.bracken_report]
        
        // report step
        software_versions = getVersions()
        workflow_params = getParams()
        report = makeReport(
            br.report,
            stats,
            software_versions.collect(), 
            workflow_params,
            template
        )

        // Kraken2 outputs with scan to accumulate results
        kr2_merged = catAssignmentsprogressive.scan(kr2.map { it -> it[1] }, taxonomy)

        // Stop server when all are processed
        stop_kraken_server(kr2_merged.assignments.collect())

        // output updating file as part of this pipeline
        output(report.report_html)
       
        output_dir(catAssignmentsprogressive.out[0].concat(progressive_kreports.out[0]))


    emit:
        report.report_html.concat(
            *outputs.collect { chan ->
                chan.flatMap { it -> [ it[1] ] }},
            software_versions,
            workflow_params
        )
}
