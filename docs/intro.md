## Introduction

wf-metagenomics offers two different approaches to assigning sequence reads to a species:

- Kraken2 offers the fastest functionality and in combination with the accompanying bracken software can be used for a more quantitative assessment of taxonomic representation within a sample.
- Minimap2 provides the finest resolution analysis but, depending on the reference database used, at the expense of significantly more compute time.

The wf-metagenomics workflow by default uses the NCBI 16S + 18S rRNA database that will be downloaded at the start of an analysis. The workflow is not tied to this database and can also be used with custom databases as required.


