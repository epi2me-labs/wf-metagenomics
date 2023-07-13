#!/usr/bin/env python
"""Parse alignments and write out a taxonomic assignment for each."""
import sys

import pandas as pd
from pysam import AlignmentFile
from .util import wf_parser  # noqa: ABS101


def main(args):
    """Write row with taxid and classification status for each alignment."""
    sam = args.sam
    output = args.output
    reference2taxid = args.ref2taxid
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
            if aln.reference_name not in ref2taxid_df.index:
                raise SystemExit(
                    """Error: The reference {} is not found in your ref2taxid file.
                    Please make sure that the ref2taxid matches the reference.
                    If your input are bam files, make sure that the ref2taxid matches
                    the reference used for the mapping step."""
                    .format(aln.reference_name))
            else:
                taxid = ref2taxid_df.at[aln.reference_name, 'taxid']

        output_tsv.write(
            '{mapped}\t{queryid}\t{taxid}\t0|{querylen}\n'.format(
                mapped=mapped, queryid=queryid, taxid=taxid, querylen=querylen
            )
        )

        aln_outfile.write(aln)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("format_minimap2")
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
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
