#!/usr/bin/env python
"""Script to run bracken or create empty output."""
import argparse
import subprocess
import sys

import pandas as pd


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument('database')
    parser.add_argument('kraken2_report')
    parser.add_argument('bracken_length')
    parser.add_argument('bracken_level')
    parser.add_argument('output')
    args = parser.parse_args()

    df_kraken2_report = pd.read_csv(
        args.kraken2_report, sep='\t',
        names=["Perc", "Count_Children", "Count_Exact", "Rank", "ID", "Name"])

    if len(df_kraken2_report.index) == 1:
        if df_kraken2_report.at[0, "Name"] == "unclassified":
            message = (
                "kraken2 report is all unclassified, "
                "writing empty output.")
            sys.stdout.write(message)
            with open(args.output, 'w') as out:
                out.write(message)
            sys.exit(0)

    ret = subprocess.run(
        'bracken '
        f'-d {args.database} '
        f'-i {args.kraken2_report} '
        f'-r {args.bracken_length} '
        f'-l {args.bracken_level} '
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


if __name__ == '__main__':
    main()
