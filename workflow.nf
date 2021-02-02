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

nextflow.enable.dsl = 2

params.help = ""
params.fastq = ""
params.out_dir = ""
params.db_path = ""

if(params.help) {
    log.info ''
    log.info 'Workflow template'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run workflow.nf [options]'
    log.info ''
    log.info 'Script Options: '
    log.info '    --fastq        FILE    Path to FASTQ file'
    log.info '    --db_path      DIR     Path centrifuge database directory'
    log.info '    --out_dir      DIR     Path for output'
    log.info ''

    return
}


process centrifuge {
    label "containerCPU"
    input:
        each file(reads)
        file db_path
    output:
        file "analysis/read_classifications.tsv"
        file "analysis/centrifuge_report.tsv"

    """
    which centrifuge
    echo "reads: $reads"
    echo "db_path: $db_path"
    centrifuge -h
    mkdir analysis
    centrifuge --met 5 --time \
        --ignore-quals -S analysis/read_classifications.tsv \
        --report-file analysis/centrifuge_report.tsv \
        -x $db_path -U $reads
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory

    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing output files"
    """
}


// workflow module
workflow pipeline {
    take:
        reads
        db_path
    main:
        test = centrifuge(reads, db_path)
    emit:
        test
}

// entrypoint workflow
workflow {
    reads = channel.fromPath(params.reads, checkIfExists:true)
    db_path = channel.fromPath(params.db_path, checkIfExists:true)
    test = pipeline(reads, db_path)
    output(test)
}
