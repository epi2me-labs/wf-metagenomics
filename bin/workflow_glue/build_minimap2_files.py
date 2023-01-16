#!/usr/bin/env python
"""Produce the ref2taxid file needed by wf-ribosomal-survey."""
import sys

import pandas as pd
from pysam import FastxFile
from .util import wf_parser  # noqa: ABS101


def main(args):
    """Write out a tsv file mapping reference names to taxids."""
    references = args.references,
    accession2taxid = args.accession2taxid,
    output_tsv = args.output_tsv,
    output_ref = args.output_ref
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


def argparser():
    """Create argument parser."""
    parser = wf_parser("build_minimap2_files")
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
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
