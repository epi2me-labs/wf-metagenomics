Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-metagenomics-report.html | Report for all samples. | aggregated |
| Abundance table with counts per taxa | abundance_table_{{ taxonomic_rank }}.tsv | Per-taxa counts TSV, including all samples. | aggregated |
| Bracken report file | bracken/{{ alias }}.kraken2_bracken.report | TSV file with the abundance of each taxa. More info about [bracken report](https://github.com/jenniferlu717/Bracken#output-kraken-style-bracken-report). | per-sample |
| Kraken2 taxonomic assignment per read (Kraken2 pipeline) | kraken2/{{ alias }}.kraken2.report.txt | Lineage-aggregated counts. More info about [kraken2 report](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format). | per-sample |
| Kraken2 taxonomic asignment per read (Kraken2 pipeline) | kraken2/{{ alias }}.kraken2.assignments.tsv | TSV file with the taxonomic assignment per read. More info about [kraken2 assignments report](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format). | per-sample |
| Host BAM file | host_bam/{{ alias }}.bam | BAM file generated from mapping filtered input reads to the host reference. | per-sample |
| BAM index file of host reads | host_bam/{{ alias }}.bai | BAM index file generated from mapping filtered input reads to the host reference. | per-sample |
| BAM file (minimap2) | bams/{{ alias }}.reference.bam | BAM file generated from mapping filtered input reads to the reference. | per-sample |
| BAM index file (minimap2) | bams/{{ alias }}.reference.bam.bai | Index file generated from mapping filtered input reads to the reference. | per-sample |
| BAM flagstat (minimap2) | bams/{{ alias }}.bamstats_results/bamstats.flagstat.tsv | Mapping results per reference | per-sample |
| Minimap2 alignment statistics (minimap2) | bams/{{ alias }}.bamstats_results/bamstats.readstats.tsv.gz | Per read stats after aligning | per-sample |
| JSON file with identified AMR genes (amr) | amr/{{ alias }}.json | JSON file with abricate results. More info about [abricate output](https://github.com/tseemann/abricate#output). | per-sample |
| Reduced reference FASTA file | igv_reference/reduced_reference.fasta.gz | Reference FASTA file containing only those sequences that have reads mapped against them. | aggregated |
| Index of the reduced reference FASTA file | igv_reference/reduced_reference.fasta.gz.fai | Index of the reference FASTA file containing only those sequences that have reads mapped against them. | aggregated |
| GZI index of the reduced reference FASTA file | igv_reference/reduced_reference.fasta.gz.gzi | Index of the reference FASTA file containing only those sequences that have reads mapped against them. | aggregated |
| JSON configuration file for IGV browser | igv.json | JSON configuration file to be loaded in IGV for visualising alignments against the reduced reference. | aggregated |
| Taxonomic assignment per read. | reads_assignments/{{ alias }}.{{kraken2|minimap2}}.assignments.tsv | TSV file with the taxonomic assignment per read. | per-sample |
