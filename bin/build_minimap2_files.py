#!/usr/bin/env python
"""Produce the ref2taxid file needed by wf-ribosomal-survey."""
import argparse
import sys

import pandas as pd
from pysam import FastxFile


def main(
    references,
    accession2taxid,
    output_tsv,
    output_ref
):
    """Write out a tsv file mapping reference names to taxids."""
    output_tsv_file = open(output_tsv, 'w')
    output_ref_file = open(output_ref, 'w')
    accession2taxid_df = pd.read_csv(
        accession2taxid, sep='\t', header=0, index_col=1)

    for ref in references:
        with FastxFile(ref) as fh:
            for entry in fh:
                try:
                    taxid = accession2taxid_df.at[entry.name, 'taxid']
                except KeyError:
                    sys.stderr.write(
                        f"Error: couldn't find taxid for {entry.name}")
                    continue
                sys.stderr.write(f"Writing {entry.name} : {taxid}\n")
                output_tsv_file.write(
                    '{name}\t{taxid}\n'.format(
                        name=entry.name, taxid=taxid))
                output_ref_file.write(str(entry) + '\n')

    output_tsv_file.close()
    output_ref_file.close()


def execute(argv):
    """Parse command line arguments and run main."""
    parser = argparse.ArgumentParser(
        description="Build a mapping of targ loci ref to taxid in .tsv format",
    )

    parser.add_argument(
        '--output_tsv',
        help=(
            "Path to output assignments in .tsv format. "
            "(Default: ref2taxid.tsv)"
        ),
        dest="output_tsv",
        default="ref2taxid.tsv",
        metavar=''
    )

    parser.add_argument(
        '--output_ref',
        help=(
            "Path to output matched sequences in .fasta format. "
            "(Default: ref2taxid.fna)"
        ),
        dest="output_ref",
        default="ref2taxid.fna",
        metavar=''
    )

    parser.add_argument(
        '-r',
        '--references',
        help="References in .fasta format.",
        dest="references",
        required=True,
        metavar='',
        nargs='*'
    )

    parser.add_argument(
        '-a',
        '--accession2taxid',
        help="Accession to taxid mapping in .tsv format.",
        dest="accession2taxid",
        required=True,
        metavar=''
    )

    args = parser.parse_args(argv)

    main(
        references=args.references,
        accession2taxid=args.accession2taxid,
        output_tsv=args.output_tsv,
        output_ref=args.output_ref
    )


if __name__ == "__main__":
    execute(sys.argv[1:])
