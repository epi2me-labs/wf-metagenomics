#!/usr/bin/env nextflow

// Developer notes
// 
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion.
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress' 

def helpMessage(){
    log.info """
wf-ribosomal-survey

Usage:
    nextflow run epi2melabs/wf-ribosomal-survey [options]

Description:
    This workflow is designed to classify reads against known organisms.

Script options:
    --fastq             DIR     Path to directory containing FASTQ files (required)
    --out_dir           DIR     Path for output (default: $params.out_dir)
    --report_name       STR     Optional report suffix (default: $params.report_name)
    --source            STR     Overrides the default reference, database and taxonomy used 
                                (Choices: ['TARGLOCI, SILVA'], Default: 'TARGLOCI')
    --taxonomy          DIR     Specifically override the taxonomy used [.tar.gz or Dir] (Default: ncbi)
    --max_len           INT     Specify read length upper limit (Default: 2000)
    --min_len           INT     Specify read length lower limit (Default: 200)

    Minimap2:
    --minimap2          BOOL    Enables classification via alignment (Default: false)
    --reference         FILE    Specifically override reference [.fna] (Default: ncbi targloci)
    --ref2taxid         FILE    Specifically override ref2taxid mapping (Default: ncbi targloci)
    --minimap2filter    STR     Filter output of minimap2 by taxids inc. child nodes, E.g. "9606,1404"
    --minimap2exclude   BOOL    Invert minimap2filter and exclude the given taxids instead
    --split_prefix      BOOL    Enable if using a very large reference with minimap2 (Default: false)

    Kraken2:
    --kraken2           BOOL    Enables classification via kmer-assignment (Default: false)
    --database          FILE    Specifically override database [.tar.gz or Dir] (Default: ncbi targeted loci)
    --kraken2filter     STR     Filter output of kraken2 by taxids inc. child nodes, E.g. "9606,1404"
    --kraken2exclude    BOOL    Invert kraken2exclude and exclude the given taxids instead
    --kraken2minimap    BOOL    Run minimap2 only on reads classified by Kraken2 (Default: true)

Notes:
    Minimap2
    The default strategy is using minimap2 to perform full
    alignments against .fasta formatted references sequences. 

    Kraken2
    It is possible to enable classification by kraken2, 
    and disable alignment which is a faster but coarser method of 
    classification reliant on the presence of a kraken2 database.

    Using both
    If both are enabled and --kraken2minimap is set, only reads 
    classified by kraken2 (and  optionally filtered by --kraken2filter) 
    are passed to the alignment step.

Notes:
    ref2taxid format is .tsv (refname  taxid), no header row.
    kraken2filter and minimap2filter accept comma sep'd taxids, e.g. 2,9606
"""
}

def prettyPrintSources(it) {
    println("> $it.key:")
    it.value.each { it2 -> {
        println(">  $it2.key:")
        it2.value.each { it3 -> {
            println(">   $it3.key: $it3.value")
        }}
    }}
}


process unpackDatabase {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        file database
    output:
        file "database_dir"
    """
    if [[ $database == *.tar.gz ]]
    then
        mkdir database_dir
        tar xf $database -C database_dir
    elif [ -d $database ]
    then
        mv $database database_dir
    else
        echo "Error: database is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}


process unpackTaxonomy {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        file taxonomy
    output:
        file "taxonomy_dir"
    """
    if [[ $taxonomy == *.tar.gz ]]
    then
        mkdir taxonomy_dir
        tar xf $taxonomy -C taxonomy_dir
    elif [ -d $taxonomy ]
    then
        mv $taxonomy taxonomy_dir
    else
        echo "Error: taxonomy is neither .tar.gz nor a dir"
        echo "Exiting".
        exit 1
    fi
    """
}


process combineFilterFastq {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        tuple path(directory), val(sample_name)
    output:
        path "${sample_name}.fastq", emit: filtered
        path "${sample_name}.stats", emit: stats
    shell:
    """
    fastcat \
        -a $params.min_len \
        -b $params.max_len \
        -q 10 \
        -s ${sample_name} \
        -r ${sample_name}.stats \
        -x ${directory} > ${sample_name}.fastq
    
    """
}


process minimap2 {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        file reads
        file reference
        file refindex
        file ref2taxid
        file taxonomy
    output:
        path "*.bam", emit: bam
        path "*.bam.bai", emit: bai
        path "*assignments.tsv", emit: assignments
        path "*lineages.txt", emit: lineage_txt
        path "*lineages.json", emit: lineage_json
    script:
        def name = reads.simpleName
        def split = params.split_prefix ? '--split-prefix tmp' : ''
    """
    minimap2 -t $params.threads -ax map-ont $split $reference $reads \
    | mapula count -r $reference -s reference -f json -p -n $name \
    | samtools view -h -F 2304 - \
    | format_minimap2.py - -o ${name}.minimap2.assignments.tsv -r $ref2taxid \
    | samtools sort -o ${name}.bam -
    samtools index ${name}.bam
    awk -F '\\t' '{print \$3}' ${name}.minimap2.assignments.tsv > taxids.tmp
    taxonkit \
        --data-dir $taxonomy \
        lineage -R taxids.tmp \
        | aggregate_lineages.py -p ${name}.minimap2
    """
}


process extractMinimap2Reads {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        file bam
        file bai
        file ref2taxid
        file taxonomy
    output:
        path "*extracted.fastq", emit: extracted
    script:
        def name = bam.simpleName
        def policy = params.minimap2exclude ? '--exclude' : ''
    """
    taxonkit \
        --data-dir $taxonomy \
        list -i $params.minimap2filter \
        --indent "" > taxids.tmp
    extract_minimap2_reads.py \
        $bam \
        -r $ref2taxid \
        -o ${name}.minimap2.extracted.fastq \
        -t taxids.tmp \
        $policy
    """
}


process kraken2 {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        file reads
        file database
        file taxonomy
    output:
        path "*.classified.fastq", emit: classified
        path "*.unclassified.fastq", emit: unclassified
        path "*lineages.txt", emit: lineage_txt
        path "*lineages.json", emit: lineage_json
        path "*assignments.tsv", emit: assignments
        path "*kraken2_report.txt", emit: kraken2_report
    script:
        def name = reads.simpleName
    """
    kraken2 \
        --db $database \
        --report ${name}.kraken2_report.txt \
        --classified-out ${name}.kraken2.classified.fastq \
        --unclassified-out ${name}.kraken2.unclassified.fastq \
        $reads > ${name}.kraken2.assignments.tsv
    awk -F '\\t' '{print \$3}' ${name}.kraken2.assignments.tsv > taxids.tmp
    taxonkit \
        --data-dir $taxonomy \
        lineage -R taxids.tmp \
        | aggregate_lineages.py -p ${name}.kraken2
    """
}


process extractKraken2Reads {
    label "wf_ribosomal_survey"
    cpus 1
    input:
        file reads
        file kraken_assignments
        file kraken_report
    output:
        path "*extracted.fastq", emit: extracted
    script:
        def taxids = (params.kraken2filter as String).replaceAll(',',' ')
        def name = reads.simpleName
        def policy = params.kraken2exclude ? '--exclude' : ''
    """
    extract_kraken_reads.py \
        -k $kraken_assignments \
        -r $kraken_report \
        -s1 $reads \
        -o ${name}.kraken2.extracted.fastq \
        -t $taxids \
        --fastq-output \
        --include-children
        $policy
    """
}


process getVersions {
    label "wf_ribosomal_survey"
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
    label "wf_ribosomal_survey"
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
    label "wf_ribosomal_survey"
    input:
        path stats
        path "versions/*"
        path "params.json"
        path lineages
    output:
        path "wf-template-*.html"
    script:
        report_name = "wf-template-" + params.report_name + '.html'
    """
    report.py \
        $report_name \
        --versions versions \
        --params params.json \
        --summaries $stats \
        --lineages $lineages
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wf_ribosomal_survey"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    take:
        samples
        reference
        refindex
        ref2taxid
        taxonomy
        database
    main:
        outputs = []
        taxonomy = unpackTaxonomy(taxonomy)
        if (params.kraken2) {
            database = unpackDatabase(database)
        }

        // Initial reads QC
        reads = combineFilterFastq(samples)

        // Run Kraken2
        if (params.kraken2) {
            kr2 = kraken2(
                reads.filtered, 
                database, 
                taxonomy
            )
            reads_to_align = kr2.classified
            outputs += [
                kr2.classified, 
                kr2.unclassified, 
                kr2.lineage_txt, 
                kr2.assignments, 
                kr2.kraken2_report
            ]
            if (params.kraken2filter) {
                kr2_filt = extractKraken2Reads(
                    reads_to_align, 
                    kr2.assignments, 
                    kr2.kraken2_report
                )
                outputs += [kr2_filt.extracted]
                reads_to_align = kr2_filt.extracted
            }
        } else {
            reads_to_align = reads.filtered
        }

        // Run Minimap2
        if (params.minimap2) {
            mm2 = minimap2(
                reads_to_align, 
                reference, 
                refindex, 
                ref2taxid, 
                taxonomy
            )
            outputs += [
                mm2.bam, 
                mm2.bai, 
                mm2.assignments, 
                mm2.lineage_txt
            ]
            if (params.minimap2filter) {
                mm2_filt = extractMinimap2Reads(
                    mm2.bam, 
                    mm2.bai, 
                    ref2taxid, 
                    taxonomy
                )
                outputs += [mm2_filt.extracted]
            }
        }

        if (params.kraken2) {
            lineages = kr2.lineage_json.collect()
        }
        if (params.minimap2) {
            lineages = mm2.lineage_json.collect()
        }

        // Reporting
        software_versions = getVersions()
        workflow_params = getParams()
        report = makeReport(
            reads.stats.collect(),
            software_versions.collect(), 
            workflow_params,
            lineages
        )
        outputs += [software_versions, workflow_params]
    emit:
        report.concat(*outputs)
}


// entrypoint workflow
workflow {
    dataDir = projectDir + '/data'

    // Ready the optional file
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    if (params.help) {
        helpMessage()
        exit 1
    }

    // Checking user parameters
    println("=================================")
    println("Checking inputs")

    if (!params.fastq) {
        helpMessage()
        println("")
        println("`--fastq` is required")
        exit 1
    }

    if (!(params.minimap2 || params.kraken2)) {
        println("")
        println("Nothing comes of nothing. Either minimap2 or kraken2 must be enabled.")
        exit 1
    }

    if (params.kraken2filter && !params.kraken2) {
        println("")
        println("Usage of kraken2filter requires `--kraken2`.")
        exit 1
    }

    source = params.source
    sources = params.sources
    if (!sources.containsKey(source)) {
        keys = sources.keySet()
        println("Default $params.defaultdb is not valid, must be one of $keys")
        exit 1
    }

    reference = file(sources[source]["reference"], type: "file")
    if (params.reference) {
        println("Checking custom reference exists")
        reference = file(params.reference, type: "file", checkIfExists:true)
    }

    refindex = file(sources[source]["refindex"], type: "file")
    if (params.reference) {
        println("Checking custom reference index exists")
        refindex = file(params.reference + '.fai', type: "file")
        if (!refindex.exists()) {
            refindex = file(OPTIONAL, type: "file")
        }
    }

    ref2taxid = file(sources[source]["ref2taxid"], type: "file")
    if (params.ref2taxid) {
        println("Checking custom ref2taxid mapping exists")
        ref2taxid = file(params.ref2taxid, type: "file", checkIfExists:true)
    }

    taxonomy = file(sources[source]["taxonomy"], type: "file")
    if (params.taxonomy) {
        println("Checking custom taxonomy mapping exists")
        ref2taxid = file(params.taxonomy, type: "dir", checkIfExists:true)
    }

    database = file(sources[source]["database"], type: "file")
    if (params.database) {
        println("Checking custom kraken2 database exists")
        database = file(params.database, type: "dir", checkIfExists:true)
    }


    // Print all params
    println("=================================")
    println("Summarising parameters")
    params.each { it -> {
        if ("$it.key" == "sources") {
            prettyPrintSources(it)
        } else {
            println("> $it.key: $it.value") 
        }
    }}

    samples = fastq_ingress(
        params.fastq, params.out_dir, params.samples, params.sanitize_fastq)

    results = pipeline(samples, reference, refindex, ref2taxid, taxonomy, database)
    output(results)
}
