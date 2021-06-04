nextflow.enable.dsl = 2

params.help = ""
params.db_path = ""
params.db_prefix = ""
params.s3_root = "https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/metagenomic_tutorial"

def helpMessage(){
    log.info """
        Workflow template

        Usage:
        nextflow run workflow.nf [options]

        Script Options:
            --db_path      DIR     Path centrifuge database directory
            --db_prefix    DIR   Name of the centrifuge database
    """
}

process download_database {
    label "containerPython"
    input:
        path db_path
    output:
        path "${params.db_prefix}.tar.gz"
    """
    mkdir ${params.db_prefix}
    wget -O ${params.db_prefix}.tar.gz "${params.s3_root}/${params.db_prefix}.tar.gz"
    """
}

process unzip_database {
    label "containerPython"
    input:
        path "${params.db_prefix}.tar.gz"
    output:
        path "${params.db_prefix}.*.cf"
    """
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
        zipped_files = download_database(db_path)
        centrifuge_files = unzip_database(zipped_files)
    emit:
        centrifuge_files
}

// entrypoint workflow
workflow {
    db_path = channel.fromPath(params.db_path, checkIfExists:true)
    db_file = new File(params.db_path, params.db_prefix + ".1.cf")
    centrifuge_db_files = download(db_path)
    output_database(centrifuge_db_files[0])
}
