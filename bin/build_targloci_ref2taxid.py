#!/usr/bin/env python
"""Produce the ref2taxid file needed by wf-ribosomal-survey."""
import argparse
import sys
from typing import List

import pandas as pd
from pysam import FastxFile


def main(
    references: List[str],
    accession2taxid: str,
    output: str
):
    """Write out a tsv file mapping reference names to taxids."""
    output_tsv = open(output, 'w')
    accession2taxid_df = pd.read_csv(
        accession2taxid, sep='\t', header=0, index_col=1)

    for ref in references:
        with FastxFile(ref) as fh:
            for entry in fh:
                try:
                    taxid = accession2taxid_df.at[entry.name, 'taxid']
                except KeyError:
                    print(
                        "Error: couldn't find taxid for {}".format(
                            entry.name
                        )
                    )
                    sys.exit(1)
                output_tsv.write(
                    '{name}\t{taxid}\n'.format(
                        name=entry.name, taxid=taxid
                    )
                )


def execute(argv) -> None:
    """Parse command line arguments and run main."""
    parser = argparse.ArgumentParser(
        description="Build a mapping of targ loci ref to taxid in .tsv format",
    )

    parser.add_argument(
        '-o',
        '--output',
        help=(
            "Path to output assignments in .tsv format. "
            "(Default: ref2taxid.targloci.tsv)"
        ),
        dest="output",
        default="ref2taxid.targloci.tsv",
        metavar='',
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
        metavar='',
    )

    args = parser.parse_args(argv)

    main(
        references=args.references,
        accession2taxid=args.accession2taxid,
        output=args.output
    )


if __name__ == "__main__":
    execute(sys.argv[1:])
