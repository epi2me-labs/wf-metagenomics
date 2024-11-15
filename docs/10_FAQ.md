If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-metagenomics/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 

+ *Which database is used by default?* - By default, the workflow uses the Standard-8 in kraken2 pipelines and the NCBI 16S + 18S rRNA database in the minimap2 workflow. It will be downloaded the first time the workflow is run and re-used in subsequent runs.

+ *Are more databases available?* - Other metagenomic databases (listed below) can be selected with the `database_set` parameter, but the workflow can also be used with a custom database if required (see [here](https://labs.epi2me.io/how-to-meta-offline/) for details).
    * 16S, 18S, ITS
        * ncbi_16s_18s and ncbi_16s_18s_28s_ITS:  Archaeal, bacterial and fungal 16S/18S and ITS data. There are two databases available using the data from [NCBI]https://www.ncbi.nlm.nih.gov/refseq/targetedloci/)
        * SILVA_138_1: The [SILVA](https://www.arb-silva.de/) database (version 138) is also available. Note that SILVA uses its own set of taxids, which do not match the NCBI taxids. We provide the respective taxdump files, but if you prefer using the NCBI ones, you can create them from the SILVA files ([NCBI](https://www.arb-silva.de/no_cache/download/archive/current/Exports/taxonomy/ncbi/)). As the SILVA database uses genus level, the last taxonomic rank at which the analysis is carried out is genus (`taxonomic_rank G`).
    * General databases (available only in kraken2 approaches)
        * Standard-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).
        * PlusPF-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core, protozoa and fungi. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).
        * PlusPFP-8: It contains references for Archaea, Bacteria, viral, plasmid, human, UniVec_Core, protozoa, fungi and plant. To use this database the memory available to the workflow must be slightly higher than size of the database index (8GB).

+ *How can I use Kraken2 indexes?* - There are different databases available [here](https://benlangmead.github.io/aws-indexes/k2).

+ *How can I use custom databases?* - If you want to run the workflow using your own Kraken2 database, you'll need to provide the database and an associated taxonomy dump. For a custom Minimap2 reference database, you'll need to provide a reference FASTA (or MMI) and an associated ref2taxid file. For a guide on how to build and use custom databases, take a look at our [article on how to run wf-metagenomics offline](https://labs.epi2me.io/how-to-meta-offline/).

+ *How can I run the workflow with less memory?* -
    When running in Kraken mode, you can set the `kraken2_memory_mapping` parameter if the available memory is smaller than the size of the database.

+ *How can I run the workflow offline?* - To run wf-metagenomics offline you can use the workflow to download the databases from the internet and prepare them for offline re-use later. If you want to use one of the databases supported out of the box by the workflow, you can run the workflow with your desired database and any input (for example, the test data). The database will be downloaded and prepared in a directory on your computer. Once the database has been prepared, it will be used automatically the next time you run the workflow without needing to be downloaded again. You can find advice on picking a suitable database in our [article on selecting databases for wf-metagenomics](https://labs.epi2me.io/metagenomic-databases/).

+ *Which databases are available for AMR?* - By default, ABRicate is set to search for AMR genes present in the [Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) database. Users can choose from a number of databases using the `amr_db` parameter. 

    |```amr_db``` | Database |
    |---------------|----------|
    |```resfinder```| [Resfinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)|
    |```ecoli_vf```| [E. coli virulence factors](https://github.com/phac-nml/ecoli_vf)|
    |```plasmidfinder```| [PlasmidFinder](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/src/master/)|
    |```card```| [Comprehensive Antibiotic Resistance Database](https://card.mcmaster.ca/)|
    |```argannot```| [ARG-ANNOT](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)|
    |```vfdb```| [Virulence factor DB](http://www.mgc.ac.cn/VFs/)|
    |```ncbi```| [NCBI AMRFinderPlus](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313047)|
    |```megares```| [MEGAres](https://www.meglab.org/megares/)|
    |```ecoh```| [E. coli AMR DB from SRST2](https://github.com/katholt/srst2/tree/master/data)|

+ *What does the `test_data` folder contain?* - This folder contains several small datasets, real and simulated, for testing the workflow with different parameters and use cases.

+ *Is it possible to run the real time approach in the cloud?* - No, the real time pipeline is not yet available to be run in the cloud.