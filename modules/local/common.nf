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