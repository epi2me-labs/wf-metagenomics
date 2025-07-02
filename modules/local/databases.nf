// Ready the optional file
OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

//minimap uses a reference.
process download_reference_ref2taxid {
    label "wfmetagenomics"
    cpus 1
    memory "2 GB"
    storeDir {params.store_dir ? "${params.store_dir}/${database_name}" : null }
    input:
        val database_name
        val reference_url
        // process ref2taxid here. It is weird that the user provides a custom reference file and it uses the default ref2taxid.
        // better to do both things together to avoid having a mismatch between both files
        // although it would be check by check_reference_ref2taxid().
        val ref2taxid_url
        val database_dir
    output:
        path "${database_dir}_db/${ref_basename}", emit: reference_file
        path "${database_dir}_db/${ref2taxid_basename}", emit: ref2taxid_file
    script:
        ref_basename = file(reference_url).Name
        ref2taxid_basename = file(ref2taxid_url).Name
    """
    wget '${reference_url}'
    wget '${ref2taxid_url}'
    mkdir ${database_dir}_db
    mv "${ref_basename}" ${database_dir}_db/
    mv "${ref2taxid_basename}" ${database_dir}_db/
    """
}

// check that each reference in the FASTA file contains its corresponding taxid in the ref2taxid file
process check_reference_ref2taxid {
    label "wfmetagenomics"
    cpus 1
    memory "2 GB"
    input:
        path reference
        path ref2taxid
    output:
        val true
    """
    # check if reference is compressed and use `zgrep` or `grep` accordingly
    grep_cmd=\$([[ $reference == *.gz ]] && echo zgrep || echo grep)

    if [[ \$(wc -l < "${ref2taxid}") -eq \$(\$grep_cmd -c '>' < "${reference}") ]]; then
        echo 'Match!'
    else
        echo "Error: The number of elements of the "${ref2taxid}" doesn't match the number of elements in the "${reference}"."
        echo "Please provide the corresponding "${ref2taxid}" for the "${reference}"."
        echo "Exiting".
        exit 1
    fi
    """
}


process unpack_download_kraken2_database {
    label "wfmetagenomics"
    cpus 1
    memory "2 GB"
    storeDir {params.store_dir ? "${params.store_dir}/${database_name}" : null }
    input:
        val database_name
        path database_local
        val database_url
        val database_dir
        val url_database_boolean
    output:
        path "${database_dir}_db", emit: database_dir
    script:
        String db_basename = database_local.name != "OPTIONAL_FILE" ? database_local.toString(): file(database_url).Name
    """   
    # Check if the folder is an url to fetch or a local path 
    if ${url_database_boolean}
    then
        wget '${database_url}'
    fi
    if [[ ${db_basename} == *.tar.gz ]]
        then
            mkdir ${database_dir}_db
            tar xf ${db_basename} -C ${database_dir}_db
    elif [[ ${db_basename} == *.zip ]]
        then
            mkdir ${db_basename}_db
            unzip  ${db_basename} -d ${database_dir}_db
    else
        echo "Error: database is neither .tar.gz , .zip"
        echo "Exiting".
        exit 1
    fi
    """
}


process determine_bracken_length {
    label "wfmetagenomics"
    cpus 1
    memory "2 GB"
    input:
        val database_name
        path database_dir
    output:
        path "${database_dir}.bracken_length.txt", emit: bracken_length_txt
    script:
        String bracken_length_output = "${database_dir}.bracken_length.txt"

    """
    if [[ -f "${database_dir}"/database${params.bracken_length}mers.kmer_distrib ]]
    then
        BRACKEN_LENGTH="${params.bracken_length}"
        echo \$BRACKEN_LENGTH > "${bracken_length_output}"
    elif ls "${database_dir}"/*.kmer_distrib
    then
        cd "${database_dir}"
        BRACKEN_LENGTH=\$(ls -v1 *.kmer_distrib | tail -1 | sed -e "s@^database@@" -e "s@mers.kmer_distrib@@")
        cd ..
        echo \$BRACKEN_LENGTH > ${bracken_length_output}
    else
        echo "Error: there is not a bracken distance file within the database folder."
        echo "This file is expected to have an extension ended in \'.kmer_distrib\'."
        echo "Exiting".
        exit 1
    fi
    """
}


// TAXONOMY
process download_unpack_taxonomy {
    label "wfmetagenomics"
    cpus 1
    memory "2 GB"
    storeDir {params.store_dir ? "${params.store_dir}/${database_name}" : null }
    input:
        val database_name
        path taxonomy_local
        val taxonomy_url //it can be a directory or a url
        val taxonomy_dir
        val url_boolean //check if it is a url or a local path
    output:
        path "${taxonomy_dir}_db" //add a tag to the name to avoid simlinks when an uncompress dir is given
    script:
        String tax_basename = taxonomy_local.name != "OPTIONAL_FILE" ? taxonomy_local.toString(): file(taxonomy_url).Name

    """
    # Check if the folder is an url to fetch or a local path 
    if ${url_boolean}
    then
        wget '${taxonomy_url}'
    fi
    # Uncompress db
    if [[ "${tax_basename}" == *.tar.gz ]]
        then
            mkdir ${taxonomy_dir}_db
            tar xf "${tax_basename}" -C ${taxonomy_dir}_db
    elif [[ "${tax_basename}" == *.zip ]]
        then
            mkdir ${taxonomy_dir}_db
            unzip  "${tax_basename}" -d ${taxonomy_dir}_db
    else
            echo "Error: taxonomy is neither .tar.gz, .zip"
            echo "Exiting".
            exit 1
    fi
    """
}

//SILVA database
process prepareSILVA {
    storeDir {params.store_dir ? "${params.store_dir}/${params.database_set}" : null }
    label "wfmetagenomics"
    cpus 2
    memory "2 GB"
    input:
        val bracken_length
    output:
        path "database/silva.fna", emit: reference
        path "taxonomy", emit: taxonomy
        path "database", emit: database
        path "seqid2taxid.map", emit: ref2taxid
        path "database/database*mers.kmer_distrib", emit: bracken_dist
    script:
    """
    kraken2-build --db ${params.database_set} --special silva
    bracken-build -d ${params.database_set} -t "${task.cpus}" -l ${bracken_length}
    # Move all the files following other default databases:
    mkdir database
    mv ${params.database_set}/*.k2d database/
    mv ${params.database_set}/library/silva.fna database/
    mv ${params.database_set}/seqid2taxid.map .
    mv ${params.database_set}/taxonomy taxonomy
    mv ${params.database_set}/database${bracken_length}mers.kmer_distrib database/
    """
}



workflow prepare_databases {
  take:
    // minimap2 and kraken2 uses different databases
    // source_data is null for custom databases
    source_data_taxonomy
    source_data_database
  main:
    log.info("Preparing databases.")
    //SILVA database is built in the moment to avoid distributing its files
    if (params.database_set == "SILVA_138_1" ){
        log.info("Note: SILVA TaxIDs do not match NCBI TaxIDs")
        log.info("Note: The database will be created from original files, which may make the wf run slower.")
        // Create all the database for both pipelines.
        if (params.classifier == "kraken2" && params.bracken_length){
            bracken_length = params.bracken_length
        } else if (params.classifier == "kraken2") {
            bracken_length = 1000
        } else {
            bracken_length = 1000 //not used
        }
        silva = prepareSILVA(bracken_length)
        reference_file = silva.reference
        ref2taxid_file = silva.ref2taxid
        taxonomy_dir = silva.taxonomy
        database_dir = silva.database
        if (params.taxonomic_rank == 'S') {
            log.info("Note: SILVA database does not provide species, genus is the lowest taxonomic rank that can be used.")
            taxonomic_rank = 'G'
        } else{
            taxonomic_rank = params.taxonomic_rank
        }
        if (params.classifier == 'kraken2') {
            bracken_length = determine_bracken_length(params.database_set, database_dir).bracken_length_txt
        }
        database_set = params.database_set
    } else{
        taxonomic_rank = params.taxonomic_rank
        // TAXONOMY: kraken2 and minimap2
        if (params.taxonomy) {
            log.info("Checking custom taxonomy mapping exists")
            // The user can provide a compress tar.gz or zip file or a directory.
            def tax = file(params.taxonomy)
            if (tax.isDirectory()) {
                // Directory: ready to use. Not stored
                taxonomy_dir = Channel.fromPath(file(params.taxonomy, type: "dir", checkIfExists:true)).first()
            } else {
                // Compress file: needs to unpack - store in store_dir
                if (tax.getScheme() != "file") { // it is an url file
                    taxonomy_local_path = OPTIONAL
                    taxonomy_remote_path = params.taxonomy
                    taxonomy_is_remote = true
                } else { // it is an file
                    taxonomy_local_path = file(params.taxonomy, type: "file", checkIfExists:true)
                    taxonomy_remote_path = '' // pass empty value in case it is a local compress file.
                    taxonomy_is_remote = false
                }             
                taxonomy_dir = download_unpack_taxonomy('custom', taxonomy_local_path, taxonomy_remote_path, tax.simpleName, taxonomy_is_remote)
            }
        } else{
            log.info("Using default taxonomy database.")
            // zip file that needs to be unpacked - store in store_dir
            taxonomy_remote_path = source_data_taxonomy.get("taxonomy", false)
            taxonomy_local_path = OPTIONAL // this is an url taxonomy db
            taxonomy_dir = download_unpack_taxonomy(params.database_set, taxonomy_local_path, taxonomy_remote_path, file(taxonomy_remote_path).simpleName, true)
        }
        // MINIMAP
        // Reference for minimap2: it can be a FASTA file or MMI index
        if (params.classifier == "minimap2") {
            database_dir = null  // not needed
            bracken_length = null // not needed
            if (params.reference) {
                // The user provides two files: FASTA/MMI and ref2taxid. - ready to use. Not stored
                log.info("Checking custom reference exists")
                def reference = file(params.reference)
                if (reference.isFile()) {
                    reference_file = Channel.fromPath(file(params.reference, type: "file", checkIfExists:true)).first()
                    // Check here that the user provides a custom ref2taxid for their custom reference.
                    // Provide a custom reference without its corresponding ref2taxid could be potentially misleading.
                    log.info("Checking custom ref2taxid mapping exists")
                    if (!params.ref2taxid) {
                        throw new Exception("Please provide the corresponding `--ref2taxid` file for your custom reference.")
                    } else {
                        ref2taxid_file = Channel.fromPath(file(params.ref2taxid, type: "file", checkIfExists:true)).first() 
                    }
                    database_set = 'custom'
                    // check that the number of elements of the reference (if it is a fasta)
                    // matches the number of elements of the ref2taxid.
                    if (! params.reference.endsWith(".mmi")) {
                        check_reference_ref2taxid(reference_file, ref2taxid_file)
                    }
                }
            } else {
                log.info("Using a default database.")
                // FASTA/MMI and ref2taxid. - download and ready to use - store in store_dir
                reference_url = source_data_database.get("reference", false)
                ref2taxid_url = source_data_database.get("ref2taxid", false)
                if (!reference_url) {
                    throw new Exception(
                        "Error: Source $database_set does not include a reference for "
                        + "use with minimap2, please choose another source, "
                        + "provide a custom reference or disable minimap2.")
                }
                reference = download_reference_ref2taxid(params.database_set, reference_url, ref2taxid_url, file(reference_url).simpleName)
                reference_file = reference.reference_file
                ref2taxid_file = reference.ref2taxid_file
                database_set = params.database_set
            }
        }
        // KRAKEN2
        // Reference for kraken2: *.k2d files
        if (params.classifier == "kraken2") {
            reference_file = null // not needed
            ref2taxid_file = null // not needed
            // The user can provide a compress tar.gz or zip file or a directory.
            if (params.database) {
                log.info("Checking custom kraken2 database exists")
                def database = file(params.database)
                if (database.isDirectory()) {
                    // Directory: ready to use. Not stored. It must include the bracken_dist
                    database_kraken2 = Channel.fromPath(file(database, type: "dir", checkIfExists:true)).first()
                    // Check it contains bracken_dist file.
                    if (!file(database.resolve('database*mers.kmer_distrib'), type: "file")) {
                    // Don't ask the user for one, because we need to copy it in the database folder and this is not stored in this case.
                    // Don't want to write any file in a customer's folder
                    // Additionally, this file is created within the kraken2 db when it is generated following bracken instructions
                    // And to run bracken, you provide the whole database with the bracken file included.
                        throw new Exception("Please make sure the corresponding `bracken_dist` file is within your custom database directory.\
                        If there are multiple bracken_dist files within the database directory, you can select one with `--bracken_length`.")
                    }
                    log.info("Using the bracken dist file within your custom database directory.")
                    database_dir = database_kraken2
                } else{
                    // Compress file: needs to unpack, include bracken_dist (separate file)- store in store_dir
                    // Make same assumtions that when input a dir as a database folder
                    // Bracken file is included in https://benlangmead.github.io/aws-indexes/k2 indexes
                    // Bracken uses the database folder assuming that the bracken file is inside
                    if (database.getScheme() != "file") { // it is an url file
                        database_local_path = OPTIONAL
                        database_remote_path = params.database
                        database_is_remote = true
                    } else { // it is an file
                        database_local_path = file(params.database, type: "file", checkIfExists:true)
                        database_remote_path = '' // pass empty value in case it is a local compress file.
                        database_is_remote = false
                    }     
                    database_kraken2 = unpack_download_kraken2_database('custom', database_local_path, database_remote_path, database.simpleName, database_is_remote)
                    database_dir = database_kraken2.database_dir
                }
                database_set = 'custom'
            } else {
                // Compress file: needs to unpack, includes bracken_dist (e.g. PlusPF-8, but also ncbi_16s_18s) - store in store_dir
                database_remote_path = source_data_database.get("database", false)
                log.info("Unpacking kraken2 indexes")
                if (!database_remote_path) {
                    throw new Exception(
                        "Error: Source $params.database_set does not include a database for "
                        + "use with kraken2, please choose another source, "
                        + "provide a custom reference or disable kraken2.")
                }
                database_local_path = OPTIONAL
                database_set = params.database_set
                database_kraken2 = unpack_download_kraken2_database(params.database_set, database_local_path, database_remote_path, file(database_remote_path).simpleName, true)
                database_dir = database_kraken2.database_dir
            }
            bracken_length = determine_bracken_length(database_set, database_dir).bracken_length_txt
        }
    }
    emit:
        taxonomy = taxonomy_dir
        ref2taxid = ref2taxid_file
        reference = reference_file
        database = database_dir
        bracken_length = bracken_length
        taxonomic_rank = taxonomic_rank
}