nextflow.enable.dsl = 2

params.help = ""
params.fastq = ""
params.out_dir = "output"
params.db_path = ""
params.db_prefix = ""

def helpMessage(){
    log.info """
        Workflow template

        Usage:
        nextflow run workflow.nf [options]

        Script Options:
            --fastq        FILE    Path to FASTQ file
            --db_path      DIR     Path centrifuge database directory
            --db_prefix      DIR     Name of the centrifuge database
            --out_dir        DIR      Name of output directory default "output"

    """
}

process centrifuge {
    label "containerCentrifuge"
    input:
        each file(fastq)
        file db_path
    output:
        file "analysis/read_classifications.tsv"
        file "analysis/centrifuge_report.tsv"
    """
    echo "fastq: fastq"
    echo "db_path: $db_path"
    echo "db_prefix": ${params.db_prefix}
    mkdir analysis
    mkdir analysis/fastq_bundles
    ls -lha $db_path*
    centrifuge --met 5 --time \
        --ignore-quals -S analysis/read_classifications.tsv \
        --report-file analysis/centrifuge_report.tsv \
        -x $db_path/${params.db_prefix} -U $fastq
    """
}

process fastcat {
    label "containerPython"
    input:
        each file(fastq)
    output:
        file "seqs.txt"
    """
    fastcat -r seqs.txt $fastq 1> /dev/null
    """
}

process generateMaster {
    label "containerPython"
    input:
        each file(fastq)
        file "analysis/read_classifications.tsv"
        file "seqs.txt"
    output:
        file "analysis/read_classification_master.tsv"
        file "wf-metagenomics-report.html"
    """
    fastcat -r seqs.txt $fastq > /dev/null
    generate_master_table.py analysis/read_classifications.tsv seqs.txt analysis --taxid 9606
    generate_report.py analysis/read_classification_master.tsv seqs.txt
    date
    """
}

process splitByMaster {
    label "containerPython"
    input:
        file fastq
        file "analysis/read_classification_master.tsv"
    output:
        path "analysis/fastq_bundles/*.fastq"
    """
    split_fastq_by_master.py $fastq analysis/read_classification_master.tsv analysis/fastq_bundles
    date
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "containerPython"
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
        fastq
        db_path
    main:
        results = centrifuge(fastq, db_path)
        seqs = fastcat(fastq)
        master = generateMaster(fastq, results[0], seqs)
        fastq = splitByMaster(fastq, master[0])
    emit:
        results[0]
        results[1]
        seqs
        master[0]
        master[1]
        fastq
}

// entrypoint workflow
workflow {
    if (params.help) {
        helpMessage()
        exit 1
    }
    if (!params.fastq | !params.db_path | !params.db_prefix) {
        helpMessage()
        println("")
        println("`--fastq`, `--db_path`, `--db_prefix` are required")
    }
    fastq = channel.fromPath(params.fastq, checkIfExists:true)
    db_path = channel.fromPath(params.db_path, checkIfExists:true)
    results = pipeline(fastq, db_path)
    output(results[0].concat(results[1], results[2], results[3], results[4], results[5]))
}
