#!/usr/bin/env python
"""Script to run bracken or create empty output."""
import pathlib
import subprocess
import sys

import pandas as pd
from .util import wf_parser  # noqa: ABS101


def main(args):
    """Run entry point."""
    df_kraken2_report = pd.read_csv(
        args.kraken2_report, sep='\t',
        names=["Perc", "Count_Children", "Count_Exact", "Rank", "ID", "Name"])

    if len(df_kraken2_report.index) == 1:
        if df_kraken2_report.at[0, "Name"] == "unclassified":
            message = (
                "kraken2 report is all unclassified, "
                "writing empty output.")
            sys.stdout.write(message)
            p = pathlib.Path(args.output)
            p.touch()
            sys.exit(0)

    ret = subprocess.run(
        'bracken '
        f'-d {args.database} '
        f'-i {args.kraken2_report} '
        f'-r {args.bracken_length} '
        f'-t {args.bracken_threshold} '
        f'-o {args.output}',
        capture_output=True,
        shell=True)

    stdout = str(ret.stdout)
    stderr = str(ret.stderr)
    sys.stdout.write(stdout)
    sys.stderr.write(stderr)

    if ret.returncode and 'no reads found' in stderr:
        with open(args.output, 'w') as out:
            out.write(stderr)
        sys.exit(0)

    sys.exit(ret.returncode)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("run_bracken")
    parser.add_argument('database')
    parser.add_argument('kraken2_report')
    parser.add_argument('bracken_length')
    parser.add_argument('bracken_level')
    parser.add_argument('bracken_threshold')
    parser.add_argument('output')
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
