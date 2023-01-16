#!/usr/bin/env python
"""Filter alignments by target reference taxid and write fastq."""
import sys

import pandas as pd
from pysam import AlignmentFile
from .util import wf_parser  # noqa: ABS101


def write_fastq_record(outfile, name, seq, qual):
    """Write record to file in fastq format."""
    outfile.write(''.join([
        '@{}\n{}\n+\n{}\n'.format(
            name,
            seq,
            # Todo: speed up
            ''.join(map(lambda x: chr(x+33), qual))
        )
    ]))


def main(args):
    """Run alignment taxonomic filtering."""
    sam = args.sam
    output = args.output
    taxids = args.taxids
    reference2taxid = args.ref2taxid
    exclude = args.exclude
    alignments = AlignmentFile(sam, "r")
    tax_ids = set(line.strip() for line in open(taxids))
    ref2taxid_df = pd.read_csv(
        reference2taxid, sep='\t', names=['acc', 'taxid'], index_col=0)

    outfile = open(output, 'w')

    for aln in alignments.fetch(until_eof=True):

        if not aln.is_unmapped:
            taxid = str(ref2taxid_df.at[aln.reference_name, 'taxid'])

            if taxid not in tax_ids:
                if exclude:
                    write_fastq_record(
                        outfile, aln.query_name,
                        aln.query_sequence, aln.query_qualities
                    )
                continue

        if not exclude:
            write_fastq_record(
                outfile, aln.query_name,
                aln.query_sequence, aln.query_qualities
            )

    outfile.close()


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("extract_minimap2_reads")
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

    parser.add_argument(
        '-t',
        help="TaxIDs to consider in line-delimited format",
        dest="taxids",
        required=True,
        metavar='',
    )

    parser.add_argument(
        '--exclude',
        help="Finds all reads NOT matching specified taxids instead",
        dest='exclude',
        required=False,
        action='store_true',
        default=False,
    )

    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
