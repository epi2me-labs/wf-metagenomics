#!/usr/bin/env python3
import argparse
from contextlib import ExitStack
from pathlib import Path
import sys

import pandas as pd
from pysam import FastxFile

parser = argparse.ArgumentParser(description="Split fastq based on label column")
parser.add_argument("fastq", type=Path, help="Path to fastq/fastq.tar.gz")
parser.add_argument(
    "master_table", type=Path, help="Path to master centrifuge table with label"
)
parser.add_argument(
    "out_dir", type=Path, help="Path to output dir for resulting fq"
)


def main(argv=None):
    args = parser.parse_args(argv)
    fq_path = args.fastq
    table = args.master_table
    out_dir = args.out_dir
    out_dir.mkdir(exist_ok=True)
    if not fq_path.is_file:
        raise FileNotFoundError(f"Fastq not found: {fq_path}")
    if not table.is_file:
        raise FileNotFoundError(f"Master table not found: {table}")

    df = pd.read_csv(table, sep=",", index_col=0)
    print(f"Table records: {len(df)}")
    labels = df["label"].unique()
    print(f"Labels found: {len(labels)}")
    with ExitStack() as stack:
        fastq_files = {
            label: stack.enter_context((out_dir / f"{label}.fastq").open("w"))
            for label in labels
        }
        for i, entry in enumerate(stack.enter_context(FastxFile(fq_path))):
            if not i % 500 and i:
                print(f"Processed: {i} records... (~ { i / len(df) * 100 :3.2f}% )")
            try:
                label = df.loc[[entry.name], "label"][0]
            except KeyError:
                print(f"Skipping: {entry.name} as not present in master table")
                continue
            fastq_files[label].write(f"{entry}")

    print("Complete!")


if __name__ == "__main__":
    sys.exit(main())
