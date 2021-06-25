nextflow.enable.dsl = 2

params.help = ""
params.fastq = ""
params.out_dir = "output"
params.db_path = ""
params.db_prefix = ""
params.threads = 8
params.download = ""
params.s3_root = "https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/metagenomic_tutorial"

def helpMessage(){
    log.info """
        Workflow template

        Usage:
        nextflow run workflow.nf [options]

        Script Options:
            --fastq        FILE    Path to FASTQ directory or path
            --db_path      DIR     Path centrifuge database directory
            --db_prefix    DIR   Name of the centrifuge database
            --out_dir      DIR   Name of output directory default: 'output'
            --threads      INT   Number of threads to run centrifuge with (default: 8)
            --s3_root      STR   S3 url location to search for prefix default: '${params.s3_root}'
    """
}

// Downloading databases workflow
process download_and_unzip_database {
    label "containerPython"
    input:
        path db_path
    output:
        path "${params.db_prefix}.*.cf"
    """
    mkdir ${params.db_prefix}
    echo "Downloading ${params.s3_root}/${params.db_prefix}.tar.gz"
    wget -O ${params.db_prefix}.tar.gz "${params.s3_root}/${params.db_prefix}.tar.gz"
    echo "Unzipping ${params.db_prefix}.tar.gz"
    tar -xzvf ${params.db_prefix}.tar.gz
    """
}

process output_database {
    label "containerPython"
    publishDir "${params.db_path}", mode: "copy", pattern: "*"
    input:
        file fname
    output:
        file fname
    """
    echo "Writing centrifuge database files"
    """
}

workflow download {
    take:
        db_path
    main:
        centrifuge_files = download_and_unzip_database(db_path)
    emit:
        centrifuge_files
}

// Running metagenomic workflow
process centrifuge {
    label "containerCentrifuge"
    input:
        file "seqs.fastq"
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
        --threads ${params.threads} \
        --ignore-quals -S analysis/read_classifications.tsv \
        --report-file analysis/centrifuge_report.tsv \
        -x $db_path/${params.db_prefix} -U seqs.fastq
    """
}

process fastcat {
    label "containerPython"
    input:
        each file(fastq)
    output:
        file "seqs.fastq"
        file "seqs.txt"
    """
    fastcat -xr seqs.txt $fastq 1> seqs.fastq
    """
}

process generateMaster {
    label "containerPython"
    input:
        file "analysis/read_classifications.tsv"
        file "seqs.txt"
    output:
        file "analysis/read_classification_master.tsv"
        file "wf-metagenomics-report.html"
    """
    generate_master_table.py analysis/read_classifications.tsv seqs.txt analysis --split "fungi:phylum:4751 bacteria:phylum:2 viruses:phylum:10239 else:superkingdom:" --human
    generate_report.py analysis/read_classification_master.tsv seqs.txt
    date
    """
}

process splitByMaster {
    label "containerPython"
    input:
        file "seqs.fastq"
        file "analysis/read_classification_master.tsv"
    output:
        path "analysis/fastq_bundles"
    """
    split_fastq_by_master.py seqs.fastq analysis/read_classification_master.tsv analysis/fastq_bundles
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
        seqs = fastcat(fastq)
        results = centrifuge(seqs[0], db_path)
        master = generateMaster(results[0], seqs[1])
        fastq = splitByMaster(seqs[0], master[0])
    emit:
        results[0]
        results[1]
        seqs[0]
        seqs[1]
        master[0]
        master[1]
        fastq
}

// entrypoint workflow
workflow {
    if (params.help) {
        helpMessage()
        exit 0
    } else if (params.download) {
        db_path = channel.fromPath(params.db_path, checkIfExists:true)
        db_file = new File(params.db_path, params.db_prefix + ".1.cf")
        centrifuge_db_files = download(db_path)
        output_database(centrifuge_db_files[0])
    } else {
        if (!params.fastq) {
            helpMessage()
            println("")
            println("`--fastq` is required")
            exit 1
        }
        if (!params.db_path) {
            helpMessage()
            println("")
            println("`--db_path` is required")
            exit 1
        }
        if (!params.db_prefix) {
            helpMessage()
            println("")
            println("`--db_prefix` is required")
            exit 1
        }
        db_file = new File(params.db_path, params.db_prefix + ".1.cf")
        if (!db_file.exists()){
            helpMessage()
            println("")
            println("Database files not present e.g. ($db_file)")
            exit 1
        }
        fastq = channel.fromPath(params.fastq, checkIfExists:true)
        db_path = channel.fromPath(params.db_path, checkIfExists:true)
        results = pipeline(fastq, db_path)
        output(results[0].concat(results[1], results[2], results[3], results[4], results[5], results[6]))
    }
}

