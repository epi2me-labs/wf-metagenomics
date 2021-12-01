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
            print("kraken2 report is all unclassified, writing empty output.")
            with open(args.output, 'w'):
                pass
            sys.exit(0)

    ret = subprocess.Popen(
        f'bracken '
        f'-d {args.database} '
        f'-i {args.kraken2_report} '
        f'-r {args.bracken_length} '
        f'-l {args.bracken_level} '
        f'-o {args.output}',
        shell=True
    ).wait()

    sys.exit(ret)


if __name__ == '__main__':
    main()
