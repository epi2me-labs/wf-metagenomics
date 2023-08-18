import groovy.json.JsonBuilder

process prepareSILVA {
    storeDir "${params.store_dir}/${params.database_set}"
    label "wfmetagenomics"
    cpus 2
    output:
        path "database/silva.fna", emit: reference
        path "taxonomy", emit: taxonomy
        path "database", emit: database
        path "seqid2taxid.map", emit: ref2taxid
        path "database1000mers.kmer_distrib", emit: bracken_dist
    script:
    """
    kraken2-build --db $params.database_set --special silva
    bracken-build -d $params.database_set -t "${task.cpus}" -l 1000
    BRACKEN_PATH=\$(dirname \$(which bracken))
    python \$BRACKEN_PATH/src/generate_kmer_distribution.py -i $params.database_set/database1000mers.kraken -o database1000mers.kmer_distrib
    # Move all the files following other default databases:
    mkdir database
    mv $params.database_set/*.k2d database/
    mv $params.database_set/library/silva.fna database/
    mv $params.database_set/seqid2taxid.map .
    mv $params.database_set/taxonomy taxonomy
    """
}

process abricateVersion {
    label "amr"
    cpus 1
    input:
        path "input_versions.txt"
    output:
        path "versions.txt"
    script:
    """
    cat input_versions.txt >> versions.txt
    abricate --version | sed 's/ /,/' >> "versions.txt"
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

workflow run_common {
    main:
        common_versions = getVersions()
        if (params.amr){
            versions = abricateVersion(common_versions)
        } else{
            versions = common_versions
        }
        parameters = getParams()
    emit:
        software_versions = versions
        parameters = parameters
}
