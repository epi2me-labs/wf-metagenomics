Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-metagenomics-report.html | Report for all samples. | aggregated |
| Abundance table with counts per taxa | ./abundance_table_{{ taxonomic_rank }}.tsv | Per-taxa counts TSV, including all samples. | aggregated |
| Bracken report file | ./bracken/{{ alias }}.kraken2_bracken.report | TSV file with the abundance of each taxa. See more info here: https://github.com/jenniferlu717/Bracken#output-kraken-style-bracken-report. | per-sample |
| Kraken2 taxonomic assignment per read (Kraken2 pipeline) | ./kraken2/{{ alias }}.kraken2.report.txt | Lineage-aggregated counts. See more info here: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format. | per-sample |
| Kraken2 taxonomic asignment per read (Kraken2 pipeline) | ./kraken2/{{ alias }}.kraken2.assignments.tsv | TSV file with the taxonomic assignment per read. See more info here: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#standard-kraken-output-format. | per-sample |
| Host BAM file | ./host_bam/{{ alias }}.bam | BAM file generated from mapping filtered input reads to the host reference. | per-sample |
| BAM index file of host reads | ./host_bam/{{ alias }}.bai | BAM index file generated from mapping filtered input reads to the host reference. | per-sample |
| BAM file (minimap2) | ./bam/{{ alias }}.bam | BAM file generated from mapping filtered input reads to the reference. | per-sample |
| BAM index file (minimap2) | ./bam/{{ alias }}.bai | Index file generated from mapping filtered input reads to the reference. | per-sample |
| JSON file with identified AMR genes (amr) | ./amr/{{ alias }}.json | JSON file with abricate results. See more info here: https://github.com/tseemann/abricate#output. | per-sample |
