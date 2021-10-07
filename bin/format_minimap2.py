#!/usr/bin/env python
"""Parse alignments and write out a taxonomic assignment for each."""
import argparse
import sys

import pandas as pd
from pysam import AlignmentFile


def main(
    sam: str,
    output: str,
    reference2taxid: str
) -> None:
    """Write row with taxid and classification status for each alignment."""
    aln_infile = AlignmentFile(sam, "r")
    aln_outfile = AlignmentFile('-', "w", template=aln_infile)
    ref2taxid_df = pd.read_csv(
        reference2taxid, sep='\t', names=['acc', 'taxid'], index_col=0)
    output_tsv = open(output, 'w+')

    for aln in aln_infile.fetch(until_eof=True):

        mapped = 'U' if aln.is_unmapped else 'C'
        queryid = aln.query_name
        querylen = aln.query_length

        taxid = 0
        if not aln.is_unmapped:
            taxid = ref2taxid_df.at[aln.reference_name, 'taxid']

        output_tsv.write(
            '{mapped}\t{queryid}\t{taxid}\t0|{querylen}\n'.format(
                mapped=mapped, queryid=queryid, taxid=taxid, querylen=querylen
            )
        )

        aln_outfile.write(aln)


def execute(argv) -> None:
    """Parse command line arguments and run main."""
    parser = argparse.ArgumentParser(
        description="Outputs assignments in a kraken2-like format",
    )

    parser.add_argument(
        help=(
            "Input alignments in SAM format. (Default: [stdin])."
        ),
        dest="sam",
        default=sys.stdin,
        metavar='',
    )

    parser.add_argument(
        '-o',
        help="Path to output assignments in .tsv format.",
        dest="output",
        required=True,
        metavar='',
    )

    parser.add_argument(
        '-r',
        help="Reference to taxid mapping in .csv format",
        dest="ref2taxid",
        required=True,
        metavar='',
    )

    args = parser.parse_args(argv)

    main(
        sam=args.sam,
        output=args.output,
        reference2taxid=args.ref2taxid,
    )


if __name__ == "__main__":
    execute(sys.argv[1:])
