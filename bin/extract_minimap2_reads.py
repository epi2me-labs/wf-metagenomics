#!/usr/bin/env python
"""Script to check that sample sheet is well-formatted."""
import sys
import argparse
import pandas as pd
from pysam import AlignmentFile


def write_fastq_record(outfile, name, seq, qual):
    outfile.write(''.join([
        '@{}\n{}\n+\n{}\n'.format(
            name,
            seq,
            # Todo: speed up
            ''.join(map(lambda x: chr(x+33), qual))
        )
    ]))


def main(
    sam: str,
    output: str,
    taxids: str,
    reference2taxid: str,
    exclude: bool = False,
) -> None:
    """
    For each alignment received, writes a formatted
    row in TSV format to the designated path which 
    contains the assigned taxid and classification
    status.
    """
    alignments = AlignmentFile(sam, "r")
    tax_ids = set(line.strip() for line in open(taxids))
    print(tax_ids)

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


def execute(argv) -> None:
    """
    Parses command line arguments and runs main.
    """
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

    args = parser.parse_args(argv)

    main(
        sam=args.sam,
        output=args.output,
        taxids=args.taxids,
        reference2taxid=args.ref2taxid,
        exclude=args.exclude,
    )


if __name__ == "__main__":
    execute(sys.argv[1:])
